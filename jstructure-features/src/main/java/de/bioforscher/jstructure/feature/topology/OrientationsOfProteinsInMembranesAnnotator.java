package de.bioforscher.jstructure.feature.topology;

import de.bioforscher.jstructure.align.AlignmentPolicy;
import de.bioforscher.jstructure.align.StructureAlignmentQuery;
import de.bioforscher.jstructure.align.impl.SingleValueDecompositionAligner;
import de.bioforscher.jstructure.model.feature.ComputationException;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.identifier.IdentifierFactory;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.nodes.Element;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Access to the internet resource of the OPM.
 * Created by bittrich on 5/17/17.
 */
@FeatureProvider(provides = { MembraneContainer.class, Topology.class })
public class OrientationsOfProteinsInMembranesAnnotator extends AbstractFeatureProvider {
    private static final Logger logger = LoggerFactory.getLogger(OrientationsOfProteinsInMembranesAnnotator.class);
    private static final String BASE_URL = "http://opm.phar.umich.edu/";
    static final String SEARCH_URL = BASE_URL + "protein.php?search=";

    @Override
    protected void processInternally(Structure protein) {
        Document document = getDocument(protein.getProteinIdentifier().getPdbId());
        processInternally(protein, document, true);
    }

    public void process(Structure protein, Document document) {
        processInternally(protein, document, false);
    }

    private void processInternally(Structure protein, Document document, boolean parseMembraneLayer) {
        try {
            MembraneContainer membraneContainer = new MembraneContainer(this);
            // extract general information - that is the first table
            Element generalDataTable = document.getElementsByClass("data").get(0);
            Element thicknessTr = generalDataTable.getElementsByTag("tr").get(1);
            membraneContainer.setHydrophobicThickness(thicknessTr.getElementsByTag("td").get(1).text());

            Element tiltTr = generalDataTable.getElementsByTag("tr").get(2);
            membraneContainer.setTiltAngle(tiltTr.getElementsByTag("td").get(1).text());

            Element transferTr = generalDataTable.getElementsByTag("tr").get(3);
            membraneContainer.setDeltaGTransfer(transferTr.getElementsByTag("td").get(1).text());

            Element topologyTr = generalDataTable.getElementsByTag("tr").get(5);
            membraneContainer.setTopology(topologyTr.getElementsByTag("td").get(1).text());

            // extract trans-membraneContainer helices - second table
            Element transMembraneSubunitsTable = document.getElementsByClass("data").get(1);
            List<TransMembraneSubunit> helices;
            try {
                helices = transMembraneSubunitsTable.getElementsByTag("tr").stream()
                        // skip the header
                        .skip(1)
                        .map(element -> element.getElementsByTag("td").get(0))
                        .map(Element::text)
                        .map(TransMembraneSubunit::new)
                        .collect(Collectors.toList());
            } catch (ArrayIndexOutOfBoundsException e) {
                // happens for pdbId: '2k21', a structure without TM-segments but rather a statement
                helices = new ArrayList<>();
            }
            membraneContainer.setTransMembraneHelices(helices);

            if (parseMembraneLayer) {
                // create global membrane object - 3rd link points to download
                String downloadLink = document.getElementById("caption").getElementsByTag("a").get(2).attr("href");
                try (InputStreamReader inputStreamReader = new InputStreamReader(new URL(BASE_URL + downloadLink).openStream())) {
                    try (BufferedReader bufferedReader = new BufferedReader(inputStreamReader)) {
                        byte[] bytes = bufferedReader.lines()
                                .collect(Collectors.joining(System.lineSeparator()))
                                .getBytes();

                        // parse protein
                        Structure opmProtein = StructureParser.source(new ByteArrayInputStream(bytes))
                                .forceProteinName(IdentifierFactory.createProteinIdentifier(downloadLink.split("=")[0].split("/")[1].substring(0, 4)))
                                .parse();

                        // superimpose opm protein onto instance of the original protein
                        //TODO this alignment is by no means perfect, but works for a first glance
                        //TODO alpha-carbon-only option
                        StructureAlignmentQuery query = StructureAlignmentQuery.of(protein, opmProtein)
                                .matchingBehavior(AlignmentPolicy.MatchingBehavior.aminoAcidsAlphaCarbonsTolerant)
                                .manipulationBehavior(AlignmentPolicy.ManipulationBehavior.COPY);
                        new SingleValueDecompositionAligner().align(query)
                                .getTransformation()
                                .transform(opmProtein);

                        // extract dummy atoms and move them to membraneContainer object
                        List<double[]> membraneAtoms = opmProtein.atoms()
                                .map(Atom::getParentGroup)
                                .filter(group -> group.getThreeLetterCode().equals("DUM"))
                                .flatMap(Group::atoms)
                                .map(Atom::getCoordinates)
                                .collect(Collectors.toList());
                        membraneContainer.setMembraneAtoms(membraneAtoms);

//                    //TODO remove, used to evaluate alignment manually
//                    Files.write(Paths.get(System.getProperty("user.home") + "/ori.pdb"), protein.getPdbRepresentation().getBytes());
//                    Files.write(Paths.get(System.getProperty("user.home") + "/opm.pdb"), opmProtein.getPdbRepresentation().getBytes());
//                    //TODO remove, used to evaluate segment positions manually
//                    Files.write(Paths.get(System.getProperty("user.home") + "/tm.pdb"), protein.select()
//                            .residueNumber(helices.stream()
//                                    .map(TransMembraneSubunit::getSegments)
//                                    .flatMap(Collection::stream)
//                                    .collect(Collectors.toList())
//                                    .toArray(new IntegerRange[0]))
//                            .asGroupContainer()
//                            .getAtomRepresentation()
//                            .getBytes());
                    }
                }
            }

            protein.getFeatureContainer().addFeature(membraneContainer);
            protein.aminoAcids()
                    .forEach(aminoAcid -> aminoAcid.getFeatureContainer().addFeature(new Topology(this,
                            membraneContainer.isTransmembraneGroup(aminoAcid))));
        } catch (Exception e) {
            logger.warn("OPM parsing for {} failed", protein.getProteinIdentifier().getFullName());
            // sometimes it is the html's fault
//            System.err.println(document.html());
            throw new ComputationException("failed to fetch or parse OPM file", e);
        }
    }

    public static Document getDocument(String pdbId) {
        try {
            Document document = getDocumentInternal(SEARCH_URL + pdbId);

            if(document.text().contains("No matches")) {
                throw new ComputationException("did not find OPM entry for " + pdbId + " - possibly it is no membrane protein");
            }

            return document;
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    static Document getDocumentInternal(String url) throws IOException {
        logger.info("fetching OPM annotation from {}", url);
        Document document = Jsoup.connect(url).get();

        // not the results directly but rather the page forwarding to representative entries
        if(document.text().contains("Representative structure(s) of this protein: ")) {
            String href = document.getElementById("body").getElementsByTag("a").first().attr("href");
            logger.info("moving to representative structure at {}", href);
            return getDocumentInternal(BASE_URL + href);
        }

        return document;
    }
}
