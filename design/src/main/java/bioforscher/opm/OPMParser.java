package bioforscher.opm;

import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.selection.Selection;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Parser for .opm files (data dumps of the OPM database by Jsoup).
 * Created by S on 29.10.2016.
 */
public class OPMParser {
    public static final String HYDROPHOBIC_THICKNESS = "HYDROPHOBIC_THICKNESS";
    public static final String TILT_ANGLE = "TILT_ANGLE";
    public static final String DELTA_G_TRANSFER = "DELTA_G_TRANSFER";
    public static final String TOPOLOGY = "TOPOLOGY";
    public static final String NUMBER_OF_TM_STRUCTURES = "NUMBER_OF_TM_STRUCTURES";
    public static final String TM_HELIX = "TM_HELIX";

    public static void parse(final Protein protein, Path path) {
        try {
            parseInternal(protein, path);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    private static void parseInternal(final Protein protein, Path path) throws IOException{
        // stores the current context - for global annotations the values always refer to the line seen last
        String context = null;
        List<TMHelix> globalHelices = new ArrayList<>();
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
                if(!Selection.on(protein)
                        .chainName(chainId)
                        .asOptionalChain()
                        .isPresent()) {
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
                    .map(p -> new Pair<>(Selection.on(protein)
                            .chainName(chainId)
                            .residueNumber(p.getLeft())
                            .asOptionalGroup(), Selection.on(protein)
                            .chainName(chainId)
                            .residueNumber(p.getRight())
                            .asOptionalGroup()))
                    //TODO replace or handle error case
                    .filter(p -> p.getLeft().isPresent() && p.getRight().isPresent())
                    .map(p -> new TMHelix(tilt, p.getLeft().get(), p.getRight().get()))
                    .collect(Collectors.toList());

                // assign tm helices to specified residues
                for (TMHelix tmhelix : helices) {
                    globalHelices.add(tmhelix);
                    tmhelix.getGroupContainer()
                            .groups()
                            .forEach(residue -> residue.setFeature(TM_HELIX, tmhelix));
                }
            }
        }

        protein.setFeature(TM_HELIX, globalHelices);
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