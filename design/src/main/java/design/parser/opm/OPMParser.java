package design.parser.opm;

import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.structure.Protein;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import static design.parser.opm.OPMParser.FeatureNames.*;

/**
 * Parser for .opm files (data dumps of the OPM database by Jsoup).
 * Created by S on 29.10.2016.
 */
public class OPMParser {
    public enum FeatureNames {
        HYDROPHOBIC_THICKNESS,
        TILT_ANGLE,
        DELTA_G_TRANSFER,
        TOPOLOGY,
        NUMBER_OF_TM_STRUCTURES,
        TM_HELIX
    }

    public static void parse(final Protein protein, Path path) {
        try {
            parseInternal(protein, path);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    private static Protein parseInternal(final Protein protein, Path path) throws IOException{
        // stores the current context - for global annotations the values always refer to the line seen last
        String context = null;
//        System.out.println(protein.getName());
        for(String line : Files.readAllLines(path)) {
            // information directly after a context line
            if(context != null && line.contains("align=\"left\"")) {
                String value = betweenTdTags(line);
                if("Hydrophobic Thickness".equals(context)) {
                    protein.setFeature(HYDROPHOBIC_THICKNESS, value);
                } else if("Tilt Angle".equals(context)) {
                    protein.setFeature(TILT_ANGLE, value);
                } else if(context.contains("G<sub>transfer</sub>")) {
                    protein.setFeature(DELTA_G_TRANSFER, value);
                } else if("Topology".equals(context)) {
                    protein.setFeature(TOPOLOGY, value);
                } else if("Number of TM Secondary Structures".equals(context)) {
                    protein.setFeature(NUMBER_OF_TM_STRUCTURES, Integer.valueOf(value));
                }
                context = null;
                continue;
            }

            // context lines always are aligned right
            if(line.contains("align=\"right\"")) {
                context = betweenBoldTags(line);
                continue;
            } else {
                context = null;
            }

            // local tm helix annotation
            if(line.contains("Tilt:")) {
                line = betweenTdTags(line);
                final String chainId = betweenBoldTags(line);
                if(!protein.findChain(chainId).isPresent()) {
                    // some OPM files refer to chains not actually present in a non-assembled PDB structure
                    continue;
                }
                final double tilt = Double.valueOf(line.split("Tilt: ")[1].split("°")[0]);
                line = line.split("Segments: ")[1];
                List<TMHelix> helices = Arrays.stream(line.split(","))
//                    .peek(System.out::println)
                    .map(s -> s.split("\\(")[1].split("\\)")[0])
                    .map(s -> s.split("-"))
                    .map(s -> new Pair<>(Integer.valueOf(s[0].trim()), Integer.valueOf(s[1].trim())))
                    .map(p -> new Pair<>(protein.findResidue(chainId, p.getFirst()), protein.findResidue(chainId, p.getSecond())))
                    //TODO replace or handle error case
                    .filter(p -> p.getFirst().isPresent() && p.getSecond().isPresent())
                    .map(p -> new TMHelix(tilt, p.getFirst().get(), p.getSecond().get()))
                    .collect(Collectors.toList());

                // assign tm helices to specified residues
                for (TMHelix tmhelix : helices) {
                    tmhelix.residues()
                           .forEach(residue -> residue.setFeature(TM_HELIX, tmhelix));
                }
            }
        }

        return protein;
    }

    /**
     * Interesting values are mostly within bold tags.
     * @param line the raw line
     * @return the first element surrounded by bold tags
     */
    private static String betweenBoldTags(String line) {
        return line.split("<b>")[1].split("</b>")[0];
    }

    /**
     * Information is written within td tags.
     * @param line the raw line
     * @return the String contained
     */
    private static String betweenTdTags(String line) {
        int begin = line.indexOf("\">") + 2;
        int end = line.indexOf("</td>");
        return line.substring(begin, end);
    }
}
