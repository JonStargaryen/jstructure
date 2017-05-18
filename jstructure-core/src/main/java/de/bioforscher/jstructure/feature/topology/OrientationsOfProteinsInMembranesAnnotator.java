package de.bioforscher.jstructure.feature.topology;

import de.bioforscher.jstructure.alignment.structure.SVDSuperimposer;
import de.bioforscher.jstructure.feature.ComputationException;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.identifier.PdbId;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
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
@FeatureProvider(provides = { GenericMembrane.class, TransMembraneHelixContainer.class })
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
                            .forceProteinName(PdbId.createFromPdbId(downloadLink.split("=")[0].split("/")[1].substring(0, 4)))
                            .parse();

                    // superimpose opm protein onto instance of the original protein
                    AtomContainer superimposedOpmProtein = SVDSuperimposer.ALPHA_CARBON_SVD_INSTANCE.align(protein, opmProtein).getAlignedQuery();

                    // extract dummy atoms and move them to membrane object
                    List<double[]> membraneAtoms = superimposedOpmProtein.atoms()
                            .map(Atom::getParentGroup)
                            .filter(group -> group.getThreeLetterCode().equals("DUM"))
                            .flatMap(Group::atoms)
                            .map(Atom::getCoordinates)
                            .collect(Collectors.toList());
                    GenericMembrane membrane = new GenericMembrane(this, membraneAtoms);
                    protein.getFeatureContainer().addFeature(membrane);
                }
            }

            // extract trans-membrane helices
            //TODO implement parsing of web page
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
