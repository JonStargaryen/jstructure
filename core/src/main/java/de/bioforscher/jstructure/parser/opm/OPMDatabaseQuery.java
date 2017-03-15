package de.bioforscher.jstructure.parser.opm;

import de.bioforscher.jstructure.feature.topology.ANVIL;
import de.bioforscher.jstructure.feature.topology.Membrane;
import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.nodes.Element;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.net.URL;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Queries the OPM database for the orientation of membrane proteins.
 * Created by bittrich on 3/13/17.
 */
public class OPMDatabaseQuery {
    private static final Logger logger = LoggerFactory.getLogger(OPMDatabaseQuery.class);
    private static final String BASE_URL = "http://opm.phar.umich.edu/";
    private static final String SEARCH_URL = BASE_URL + "protein.php?search=";
    public static final String HOMOLOGOUS_PROTEINS = "HOMOLOGOUS_PROTEINS";

    public static Protein parseAnnotatedProteinById(String pdbId) {
        try {
            Document document = getDocument(SEARCH_URL + pdbId);
            // 3rd link points to download
            String downloadLink = document.getElementById("caption").getElementsByTag("a").get(2).attr("href");

            try(InputStreamReader inputStreamReader = new InputStreamReader(new URL(BASE_URL + downloadLink).openStream())) {
                try (BufferedReader bufferedReader = new BufferedReader(inputStreamReader)) {
                    byte[] bytes = bufferedReader.lines()
                            .collect(Collectors.joining(System.lineSeparator()))
                            .getBytes();

                    // parse protein
                    Protein protein = ProteinParser.source(new ByteArrayInputStream(bytes))
                            .forceProteinName(downloadLink.split("=")[0].split("/")[1].substring(0, 4))
                            .parse();

                    // retrieve homologous protein ids (if any)
                    List<String> homologous = document.getElementsByClass("data")
                            .first()
                            .getElementsByClass("row1")
                            .get(3)
                            .getElementsByTag("a")
                            .stream()
                            .map(Element::text)
                            .distinct()
                            .collect(Collectors.toList());
                    protein.setFeature(HOMOLOGOUS_PROTEINS, homologous);

                    // parse dummy membrane proteins and assign them in standardized format
                    List<Group> membraneGroups = protein.chains()
                            .flatMap(Chain::groups)
                            .filter(group -> group.getThreeLetterCode().equals("DUM"))
                            .collect(Collectors.toList());
                    protein.chains().forEach(chain -> chain.getGroups().removeIf(membraneGroups::contains));

                    Membrane membrane = new Membrane(membraneGroups.stream()
                            .map(Group::getAtoms)
                            .flatMap(Collection::stream)
                            .map(Atom::getCoordinates)
                            .filter(atom -> minimalSquaredDistanceToProteinAtom(protein, atom) > 12.0)
                            .collect(Collectors.toList()));
                    protein.setFeature(ANVIL.MEMBRANE, membrane);

                    return protein;
                }
            }
        } catch (IOException e) {
            throw new UncheckedIOException("failed to fetch OPM files", e);
        }
    }

    /**
     * Computes the minimal distance of a point to any alpha carbon of a structure.
     * @param membranePseudoAtom the point to check
     * @return the closest distance to any alpha carbon of the protein
     */
    private static double minimalSquaredDistanceToProteinAtom(Protein protein, final double[] membranePseudoAtom) {
        return protein.atoms()
        .map(Atom::getCoordinates)
        .mapToDouble(coordinates -> LinearAlgebra3D.distanceFast(coordinates, membranePseudoAtom))
        .min()
        .orElse(Double.MAX_VALUE);
    }
    
    static Document getDocument(String url) throws IOException {
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
