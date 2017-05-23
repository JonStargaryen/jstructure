package de.bioforscher.jstructure.feature.topology;

import de.bioforscher.jstructure.alignment.SVDSuperimposer;
import de.bioforscher.jstructure.feature.ComputationException;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.identifier.ProteinIdentifier;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.nodes.Element;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Access to the internet resource of the OPM.
 * Created by bittrich on 5/17/17.
 */
@FeatureProvider(provides = Membrane.class)
public class OrientationsOfProteinsInMembranesAnnotator extends AbstractFeatureProvider {
    private static final Logger logger = LoggerFactory.getLogger(OrientationsOfProteinsInMembranesAnnotator.class);
    private static final String BASE_URL = "http://opm.phar.umich.edu/";
    static final String SEARCH_URL = BASE_URL + "protein.php?search=";

    @Override
    protected void processInternally(Protein protein) {
        try {
            Document document = getDocument(SEARCH_URL + protein.getPdbId().getPdbId());

            if(document.text().contains("No matches")) {
                throw new ComputationException("did not find OPM entry for " + protein.getIdentifier() + " - possibly it is no membrane protein");
            }

            // create global membrane object - 3rd link points to download
            String downloadLink = document.getElementById("caption").getElementsByTag("a").get(2).attr("href");
            try (InputStreamReader inputStreamReader = new InputStreamReader(new URL(BASE_URL + downloadLink).openStream())) {
                try (BufferedReader bufferedReader = new BufferedReader(inputStreamReader)) {
                    byte[] bytes = bufferedReader.lines()
                            .collect(Collectors.joining(System.lineSeparator()))
                            .getBytes();

                    // parse protein
                    Protein opmProtein = ProteinParser.source(new ByteArrayInputStream(bytes))
                            .forceProteinName(ProteinIdentifier.createFromPdbId(downloadLink.split("=")[0].split("/")[1].substring(0, 4)))
                            .parse();
                    Membrane membrane = new Membrane(this);

                    // superimpose opm protein onto instance of the original protein
                    //TODO this alignment is by no means perfect, but works for a first glance
                    SVDSuperimposer.ALPHA_CARBON_SVD_INSTANCE
                            .align(protein.select()
                                    .aminoAcids()
                                    .asGroupContainer(),
                                   opmProtein.select()
                                    .aminoAcids()
                                    .asGroupContainer())
                            .transform(opmProtein);

                    // extract dummy atoms and move them to membrane object
                    List<double[]> membraneAtoms = opmProtein.atoms()
                            .map(Atom::getParentGroup)
                            .filter(group -> group.getThreeLetterCode().equals("DUM"))
                            .flatMap(Group::atoms)
                            .map(Atom::getCoordinates)
                            .collect(Collectors.toList());
                    membrane.setMembraneAtoms(membraneAtoms);

                    // extract general information - that is the first table
                    Element generalDataTable = document.getElementsByClass("data").get(0);
                    Element thicknessTr = generalDataTable.getElementsByTag("tr").get(1);
                    membrane.setHydrophobicThickness(thicknessTr.getElementsByTag("td").get(1).text());

                    Element tiltTr = generalDataTable.getElementsByTag("tr").get(2);
                    membrane.setTiltAngle(tiltTr.getElementsByTag("td").get(1).text());

                    Element transferTr = generalDataTable.getElementsByTag("tr").get(3);
                    membrane.setDeltaGTransfer(transferTr.getElementsByTag("td").get(1).text());

                    Element topologyTr = generalDataTable.getElementsByTag("tr").get(5);
                    membrane.setTopology(topologyTr.getElementsByTag("td").get(1).text());

                    // extract trans-membrane helices - second table
                    Element transMembraneSubunitsTable = document.getElementsByClass("data").get(1);
                    List<TransMembraneHelix> helices = transMembraneSubunitsTable.getElementsByTag("tr").stream()
                            // skip the header
                            .skip(1)
                            .map(element -> element.getElementsByTag("td").get(0))
                            .map(Element::text)
                            .map(TransMembraneHelix::new)
                            .collect(Collectors.toList());
                    membrane.setTransMembraneHelices(helices);

                    protein.getFeatureContainer().addFeature(membrane);
                    //TODO remove, used to evaluate alignment manually
//                    Files.write(Paths.get(System.getProperty("user.home") + "/ori.pdb"), protein.getPdbRepresentation().getBytes());
//                    Files.write(Paths.get(System.getProperty("user.home") + "/opm.pdb"), opmProtein.getPdbRepresentation().getBytes());
                }
            }
        } catch (IOException e) {
            throw new ComputationException("failed to fetch OPM file", e);
        }
    }

    Document getDocument(String url) throws IOException {
        logger.info("fetching OPM annotation from {}", url);
        Document document = Jsoup.connect(url).get();

        // not the results directly but rather the page forwarding to representative entries
        if(document.text().contains("Representative structure(s) of this protein: ")) {
            String href = document.getElementById("body").getElementsByTag("a").first().attr("href");
            logger.info("moving to representative structure at {}", href);
            return getDocument(BASE_URL + href);
        }

        return document;
    }
}
