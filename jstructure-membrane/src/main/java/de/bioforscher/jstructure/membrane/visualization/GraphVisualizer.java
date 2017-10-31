package de.bioforscher.jstructure.membrane.visualization;

import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.ObjectMapper;
import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceAreaCalculator;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfileCalculator;
import de.bioforscher.jstructure.feature.sse.dssp.DictionaryOfProteinSecondaryStructure;
import de.bioforscher.jstructure.mathematics.graph.PartitionedGraph;
import de.bioforscher.jstructure.mathematics.graph.partitioning.Module;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.util.*;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class GraphVisualizer {
    // http://tools.medialab.sciences-po.fr/iwanthue/ - 20 distinct colors
    private static final List<String> PALETTE = Pattern.compile(", ").splitAsStream("[207,61,125], " +
            "[94,198,68], " +
            "[107,61,202], " +
            "[106,156,55], " +
            "[202,74,206], " +
            "[183,171,53], " +
            "[105,111,209], " +
            "[213,136,66], " +
            "[112,46,133], " +
            "[88,181,125], " +
            "[212,73,57], " +
            "[81,169,167], " +
            "[125,47,54], " +
            "[116,148,198], " +
            "[126,85,43], " +
            "[203,122,194], " +
            "[63,95,47], " +
            "[87,64,110], " +
            "[171,165,109], " +
            "[206,134,145]").collect(Collectors.toList());

    public static String composeJsonRepresentation(Chain chain,
                                                   List<AminoAcid> earlyFoldingResidues,
                                                   List<Double> dynamineScores,
                                                   List<Double> efoldmineScores,
                                                   PartitionedGraph<AminoAcid> experimentalData,
                                                   Map<String, PartitionedGraph<AminoAcid>> inSilicoData) throws JsonProcessingException {
        // annotate secondary structure elements
        new DictionaryOfProteinSecondaryStructure().process(chain.getParentStructure());
        new EnergyProfileCalculator().process(chain.getParentStructure());
        new AccessibleSurfaceAreaCalculator().process(chain.getParentStructure());

        return new ObjectMapper().writeValueAsString(new JsonChain(chain,
                earlyFoldingResidues,
                dynamineScores,
                efoldmineScores,
                experimentalData,
                inSilicoData));
    }

    public static String composeDefinitionString(PartitionedGraph<AminoAcid> graph) {
        Map<Optional<Module<AminoAcid>>, List<AminoAcid>> map = graph.getNodes().stream()
                .collect(Collectors.groupingBy(graph::getModuleOf));

        return map.entrySet().stream()
                .filter(entry -> entry.getKey().isPresent())
                .map(GraphVisualizer::composeDefinitionString)
                .collect(Collectors.joining(System.lineSeparator()));
    }

    private static String composeDefinitionString(Map.Entry<Optional<Module<AminoAcid>>, List<AminoAcid>> entry) {
        Module<AminoAcid> module = entry.getKey().orElseThrow(() -> new NoSuchElementException(""));
        String id = module.getIdentifier();

        List<AminoAcid> aminoAcids = entry.getValue().stream()
                .map(AminoAcid.class::cast)
                .collect(Collectors.toList());

        StringJoiner ranges = new StringJoiner(",");
        String currentRange = "";
        String lastId = "";
        for(AminoAcid aminoAcid : aminoAcids) {
            String currentId = aminoAcid.getResidueIdentifier().toString();
            Optional<AminoAcid> previousAminoAcid = aminoAcid.getPreviousAminoAcid();
            if(!previousAminoAcid.isPresent() || !module.containsNode(previousAminoAcid.get())) {
                if (!currentRange.equals("")) {
                    // terminate last range record if something was observed before
                    currentRange = currentRange + "-" + lastId;
                    ranges.add(currentRange);
                }
                currentRange = currentId;
            }

            lastId = currentId;
        }
        ranges.add(currentRange + "-" + lastId);

        return id + ":" + ranges.toString();
    }

    public static String getPyMolString(PartitionedGraph<AminoAcid> graph) {
        AminoAcid aminoAcid = graph.getNodes().get(0);
        String pdbId = aminoAcid.getParentChain().getChainIdentifier().getProteinIdentifier().getPdbId();
        String chainId = aminoAcid.getParentChain().getChainIdentifier().getChainId();

        return graph.getNodes().stream()
                .map(node -> composeColorString(graph, node))
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.joining(System.lineSeparator(),
                        "fetch " + pdbId + ", async=0" + System.lineSeparator() +
                                "bg_color white" + System.lineSeparator() +
                                // hide non-relevant stuff
                                "hide everything" + System.lineSeparator() +
                                "show cartoon, chain " + chainId + System.lineSeparator() +
                                // decolor everything
                                "color grey80" + System.lineSeparator() +
                                // create custom color palette
                                IntStream.range(0, PALETTE.size())
                                        .mapToObj(i -> {
                                            String[] split = PALETTE.get(i)
                                                    .replace("[", "")
                                                    .replaceAll("]", "")
                                                    .split(",");
                                            return "set_color cc" + i + "=[" +
                                                    Double.valueOf(split[0]) / 255 + ", " +
                                                    Double.valueOf(split[1]) / 255 + ", " +
                                                    Double.valueOf(split[2]) / 255 + "]";
                                        })
                                        .collect(Collectors.joining(System.lineSeparator())) + System.lineSeparator(),
                        ""));
    }

    private static Optional<String> composeColorString(PartitionedGraph<AminoAcid> graph, AminoAcid node) {
        Optional<Module<AminoAcid>> module = graph.getModuleOf(node);

        // ignore unassigned residues
        if(!module.isPresent()) {
            return Optional.empty();
        }

        // select and set color
        int moduleNumber = graph.getModules().indexOf(module.get()) % PALETTE.size();
        //TODO insertion codes
        int resi = node.getResidueIdentifier().getResidueNumber();

        return Optional.of("color cc" + moduleNumber + ", resi " + resi);
    }
}
