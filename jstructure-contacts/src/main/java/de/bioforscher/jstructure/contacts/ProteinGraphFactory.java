package de.bioforscher.jstructure.contacts;

import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.sse.dssp.DSSPSecondaryStructure;
import de.bioforscher.jstructure.feature.sse.dssp.DictionaryOfProteinSecondaryStructure;
import de.bioforscher.jstructure.mathematics.graph.Edge;
import de.bioforscher.jstructure.mathematics.graph.Graph;
import de.bioforscher.jstructure.mathematics.graph.PartitionedGraph;
import de.bioforscher.jstructure.mathematics.graph.partitioning.Module;
import de.bioforscher.jstructure.mathematics.graph.partitioning.impl.MarkovClusteringAlgorithm;
import de.bioforscher.jstructure.mathematics.SetOperations;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.function.BiPredicate;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Constructs graph representations of protein chains.
 */
public class ProteinGraphFactory {
    private static final double COVALENT_BOND_WEIGHT = 10;
    private static final double SELF_LOOP_WEIGHT = 10;
    private static final double CONTACT_WEIGHT = 5;

//    public enum WeightingScheme {
//        UNWEIGHTED(interaction -> 1.0),
//        CLASSIFIED(interaction -> {
//            if(interaction instanceof HydrogenBond) {
//                return 5;
//            } else if(interaction instanceof HydrophobicInteraction) {
//                return 2;
//            } else {
//                // covalent bond or any of the freak interactions
//                return 10;
//            }
//        }),
//        CLASSIFIED(interaction -> {
//            if(interaction instanceof HalogenBond) {
//                return 2;
//            } else if(interaction instanceof HydrogenBond) {
//                return 3;
//            } else if(interaction instanceof HydrophobicInteraction) {
//                return 1;
//            } else if(interaction instanceof MetalComplex) {
//                return 3;
//            } else if(interaction instanceof PiCationInteraction) {
//                return 2;
//            } else if(interaction instanceof PiStacking) {
//                return 2;
//            } else if(interaction instanceof SaltBridge) {
//                return 2;
//            } else if(interaction instanceof WaterBridge) {
//                return 3;
//            } else {
//                // covalent bond
//                return 4;
//            }
//        }),
//        ENERGETIC(interaction -> {
//            if(interaction instanceof HalogenBond) {
//                return 8;
//            } else if(interaction instanceof HydrogenBond) {
//                return 25;
//            } else if(interaction instanceof HydrophobicInteraction) {
//                // see: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3086625/
//                return 1.1 * 4.184;
//            } else if(interaction instanceof MetalComplex) {
//                //TODO validate, find reference
//                return 7 *  4.184;
//            } else if(interaction instanceof PiCationInteraction) {
//                return 13;
//            } else if(interaction instanceof PiStacking) {
//                return ((PiStacking) interaction).getType().equals("T") ? 11 : 8;
//            } else if(interaction instanceof SaltBridge) {
//                return 8;
//            } else if(interaction instanceof WaterBridge) {
//                // see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2144129/pdf/10548043.pdf
//                return 7.4 * 4.184;
//            } else {
//                // covalent bond
//                return 40;
//            }
//        });
//
//        private final ToDoubleFunction<PLIPInteraction> weightingFunction;
//
//        WeightingScheme(ToDoubleFunction<PLIPInteraction> weightingFunction) {
//            this.weightingFunction = weightingFunction;
//        }
//
//        public double getWeight(PLIPInteraction plipInteraction) {
//            return weightingFunction.applyAsDouble(plipInteraction);
//        }
//
//        public double getCovalentWeight() {
//            return getWeight(null);
//        }
//    }

    public enum InteractionScheme {
        KRISHNAN2007((aminoAcid1, aminoAcid2) -> SetOperations.cartesianProductOf(aminoAcid1.getAtoms(), aminoAcid2.getAtoms())
                .anyMatch(pair -> pair.getLeft().calculate().distance(pair.getRight()) <
                        pair.getLeft().getElement().getVDWRadius() + pair.getRight().getElement().getVDWRadius() + 1.0)),
        KHAN2015((aminoAcid1, aminoAcid2) -> SetOperations.cartesianProductOf(aminoAcid1.getAtoms(), aminoAcid2.getAtoms())
                .anyMatch(pair -> pair.getLeft().calculate().distance(pair.getRight()) <
                        pair.getLeft().getElement().getVDWRadius() + pair.getRight().getElement().getVDWRadius() + 0.6)),
        HLEAP2013((aminoAcid1, aminoAcid2) -> SetOperations.cartesianProductOf(aminoAcid1.getAtoms(), aminoAcid2.getAtoms())
                .anyMatch(pair -> pair.getLeft().calculate().distanceFast(pair.getRight()) < 4.5 * 4.5)),
        FISCHER2000((aminoAcid1, aminoAcid2) -> SetOperations.cartesianProductOf(aminoAcid1.getAtoms(), aminoAcid2.getAtoms())
                .anyMatch(pair -> pair.getLeft().calculate().distanceFast(pair.getRight()) < 4 * 4)),
        CALPHA7((aminoAcid1, aminoAcid2) -> aminoAcid1.getCa().calculate()
                .distanceFast(aminoAcid2.getCa()) < 7 * 7),
        CALPHA8((aminoAcid1, aminoAcid2) -> aminoAcid1.getCa().calculate()
                .distanceFast(aminoAcid2.getCa()) < 8 * 8),
        CALPHA10((aminoAcid1, aminoAcid2) -> aminoAcid1.getCa().calculate()
                .distanceFast(aminoAcid2.getCa()) < 10 * 10),
        CALPHA12((aminoAcid1, aminoAcid2) -> aminoAcid1.getCa().calculate()
                .distanceFast(aminoAcid2.getCa()) < 12 * 12),
        SALENTIN2015((aminoAcid1, aminoAcid2) -> aminoAcid1.getParentChain().getFeature(PLIPInteractionContainer.class)
                .areInContact(aminoAcid1, aminoAcid2));

        private final BiPredicate<AminoAcid, AminoAcid> criterion;

        InteractionScheme(BiPredicate<AminoAcid, AminoAcid> criterion) {
            this.criterion = criterion;
        }

        public boolean areInContact(AminoAcid aminoAcid1, AminoAcid aminoAcid2) {
            return criterion.test(aminoAcid1, aminoAcid2);
        }
    }

    public static Graph<AminoAcid> createProteinGraph(Chain chain, InteractionScheme interactionScheme) {
        // ensure PLIP data is present when needed
        if (interactionScheme == InteractionScheme.SALENTIN2015) {
            if(chain.getFeatureContainer().getFeatures().stream().noneMatch(PLIPInteractionContainer.class::isInstance)) {
                throw new IllegalArgumentException("no plip intra-molecular interactions annotated - cannot graph representation");
            }
        }

        List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());
        List<Edge<AminoAcid>> interactions = new ArrayList<>();

        for(int i = 0; i < aminoAcids.size() - 1; i++) {
            AminoAcid aminoAcid1 = aminoAcids.get(i);
            for(int j = i + 1; j < aminoAcids.size(); j++) {
                AminoAcid aminoAcid2 = aminoAcids.get(j);
                if(j == i + 1) {
                    // comment next line to ignore consecutive amino acids
                    interactions.add(new Edge<>(aminoAcid1, aminoAcid2, COVALENT_BOND_WEIGHT));
                    continue;
                }

                if(interactionScheme.areInContact(aminoAcid1, aminoAcid2)) {
                    interactions.add(new Edge<>(aminoAcid1, aminoAcid2, CONTACT_WEIGHT));
                    //TODO reimpl weighting
                }
            }
        }

        return new Graph<>(aminoAcids, interactions);
    }

    public static class Partitioned {
        public static PartitionedGraph<AminoAcid> fromNetCartoFile(Graph<AminoAcid> graph, InputStream netCartoFile) {
            List<Module<AminoAcid>> modules = new BufferedReader(new InputStreamReader(netCartoFile))
                    .lines()
                    .filter(line -> !line.startsWith("#"))
                    .map(line -> new Module<>(line.split(" ")[0],
                            graph,
                            Pattern.compile("\\s+").splitAsStream(line.split("---")[1].trim())
                            .mapToInt(Integer::valueOf)
                            .mapToObj(residueNumber -> graph.getNodes()
                                    .stream()
                                    .filter(node -> node.getResidueIdentifier().getResidueNumber() == residueNumber)
                                    .findFirst())
                            .filter(Optional::isPresent)
                            .map(Optional::get)
                            .collect(Collectors.toList())))
                    .collect(Collectors.toList());
            return new PartitionedGraph<>(graph, modules);
        }

        public static PartitionedGraph<AminoAcid> fromNetCartoFile(Graph<AminoAcid> graph, Path netCartoFile) {
            try {
                return fromNetCartoFile(graph, Files.newInputStream(netCartoFile));
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
        }

        public static PartitionedGraph<AminoAcid> fromSecondaryStructureElements(Graph<AminoAcid> graph) {
            Chain chain = graph.getNodes().get(0).getParentChain();
            Structure structure = chain.getParentStructure();

            new DictionaryOfProteinSecondaryStructure().process(structure);

            List<Module<AminoAcid>> modules = new ArrayList<>();
            List<AminoAcid> currentModule = new ArrayList<>();
            String lastSse = null;

            for(AminoAcid aminoAcid : graph.getNodes()) {
                String sse = aminoAcid.getFeature(DSSPSecondaryStructure.class)
                        .getSecondaryStructure()
                        .getReducedRepresentation();
                // handle initial sse
                if(lastSse == null) {
                    lastSse = sse;
                }

                // we define modules as change of sse
                if(sse.equals(lastSse)) {
                    // when it matches, append last module
                    currentModule.add(aminoAcid);
                } else {
                    // if not: change sse, assign module and clear list
                    lastSse = sse;
                    if(!currentModule.isEmpty()) {
                        modules.add(new Module<>(String.valueOf(modules.size() + 1),
                                graph,
                                currentModule));
                    }
                    currentModule = new ArrayList<>();
                    currentModule.add(aminoAcid);
                }
            }

            modules.add(new Module<>(String.valueOf(modules.size() + 1),
                    graph,
                    currentModule));

            return new PartitionedGraph<>(graph, modules);
        }

        public static PartitionedGraph<AminoAcid> fromDefinitionFile(Graph<AminoAcid> graph, InputStream definitionFile) {
            List<Module<AminoAcid>> modules = new BufferedReader(new InputStreamReader(definitionFile))
                    .lines()
                    .filter(line -> !line.startsWith("#"))
                    .map(line -> {
                        String id = line.split(":")[0];
                        String rawRanges = line.split(":")[1];
                        List<AminoAcid> nodes = Pattern.compile(",").splitAsStream(rawRanges)
                                .map(rawRange -> rawRange.split("-"))
                                .flatMap(rawRange -> IntStream.range(Integer.valueOf(rawRange[0]), Integer.valueOf(rawRange[1]) + 1).boxed())
                                .map(residueNumber -> graph.getNodes().stream().filter(node -> node.getResidueIdentifier().getResidueNumber() == residueNumber).findFirst().get())
                                .collect(Collectors.toList());
                        return new Module<>(id,
                                graph,
                                nodes);
                    })
                    .collect(Collectors.toList());

            return new PartitionedGraph<>(graph, modules);
        }

        public static PartitionedGraph<AminoAcid> fromDefinitionFile(Graph<AminoAcid> graph, Path definitionFile) {
            try {
                return fromDefinitionFile(graph, Files.newInputStream(definitionFile));
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
        }

        public static PartitionedGraph<AminoAcid> fromMCL(Graph<AminoAcid> graph,
                                                          double expand,
                                                          double inflation) {
            return new MarkovClusteringAlgorithm(expand,
                    inflation,
                    SELF_LOOP_WEIGHT,
                    MarkovClusteringAlgorithm.DEFAULT_MAX_ITERATIONS,
                    MarkovClusteringAlgorithm.DEFAULT_SIMILARITY_FUNCTION).partitionGraph(graph);
        }
    }
}
