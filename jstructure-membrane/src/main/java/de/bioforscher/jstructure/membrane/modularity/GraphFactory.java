package de.bioforscher.jstructure.membrane.modularity;

import de.bioforscher.jstructure.feature.interactions.*;
import de.bioforscher.jstructure.feature.sse.dssp.DSSPSecondaryStructure;
import de.bioforscher.jstructure.feature.sse.dssp.DictionaryOfProteinSecondaryStructure;
import de.bioforscher.jstructure.mathematics.graph.Edge;
import de.bioforscher.jstructure.mathematics.graph.Graph;
import de.bioforscher.jstructure.mathematics.graph.PartitionedGraph;
import de.bioforscher.jstructure.mathematics.graph.partitioning.Module;
import de.bioforscher.jstructure.mathematics.graph.partitioning.algorithms.MCL;
import de.bioforscher.jstructure.model.SetOperations;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.jsoup.Jsoup;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.function.ToDoubleFunction;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Constructs graph representations of protein chains.
 */
public class GraphFactory {
    private static final Logger logger = LoggerFactory.getLogger(GraphFactory.class);
    private static final PLIPIntraMolecularAnnotator plipIntraMolecularAnnotator = new PLIPIntraMolecularAnnotator();

    public enum WeightingScheme {
        UNWEIGHTED(interaction -> 1.0),
        CLASSIFIED(interaction -> {
            if(interaction instanceof HalogenBond) {
                return 2;
            } else if(interaction instanceof HydrogenBond) {
                return 3;
            } else if(interaction instanceof HydrophobicInteraction) {
                return 1;
            } else if(interaction instanceof MetalComplex) {
                return 3;
            } else if(interaction instanceof PiCationInteraction) {
                return 2;
            } else if(interaction instanceof PiStacking) {
                return 2;
            } else if(interaction instanceof SaltBridge) {
                return 2;
            } else if(interaction instanceof WaterBridge) {
                return 3;
            } else {
                // covalent bond
                return 4;
            }
        }),
        ENERGETIC(interaction -> {
            if(interaction instanceof HalogenBond) {
                return 8;
            } else if(interaction instanceof HydrogenBond) {
                return 25;
            } else if(interaction instanceof HydrophobicInteraction) {
                // see: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3086625/
                return 1.1 * 4.184;
            } else if(interaction instanceof MetalComplex) {
                //TODO validate, find reference
                return 7 *  4.184;
            } else if(interaction instanceof PiCationInteraction) {
                return 13;
            } else if(interaction instanceof PiStacking) {
                return ((PiStacking) interaction).getType().equals("T") ? 11 : 8;
            } else if(interaction instanceof SaltBridge) {
                return 8;
            } else if(interaction instanceof WaterBridge) {
                // see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2144129/pdf/10548043.pdf
                return 7.4 * 4.184;
            } else {
                // covalent bond
                return 40;
            }
        });

        private final ToDoubleFunction<PLIPInteraction> weightingFunction;

        WeightingScheme(ToDoubleFunction<PLIPInteraction> weightingFunction) {
            this.weightingFunction = weightingFunction;
        }

        public double getWeight(PLIPInteraction plipInteraction) {
            return weightingFunction.applyAsDouble(plipInteraction);
        }

        public double getCovalentWeight() {
            return getWeight(null);
        }
    }

    public static Graph<AminoAcid> createGraphFromPlipDocument(Chain chain, String plipDocument, WeightingScheme scheme) {
        plipIntraMolecularAnnotator.process(chain, Jsoup.parse(plipDocument));

        List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());
        List<Edge<AminoAcid>> interactions = new ArrayList<>();

        for(int i = 0; i < aminoAcids.size() - 1; i++) {
            AminoAcid aminoAcid1 = aminoAcids.get(i);
            PLIPInteractionContainer interactions1 = aminoAcid1.getFeature(PLIPInteractionContainer.class);

            for(int j = i + 1; j < aminoAcids.size(); j++) {
                AminoAcid aminoAcid2 = aminoAcids.get(j);
                if(j == i + 1) {
                    // comment next line to ignore consecutive amino acids
                    interactions.add(new Edge<>(aminoAcid1, aminoAcid2, scheme.getCovalentWeight()));
                    continue;
                }

                PLIPInteractionContainer interactions2 = aminoAcid2.getFeature(PLIPInteractionContainer.class);

                Optional<PLIPInteraction> potentialInteraction = interactions1.getInteractions()
                        .stream()
                        .filter(plipInteraction -> interactions2.getInteractions().contains(plipInteraction))
                        .findFirst();
                potentialInteraction.ifPresent(plipInteraction ->
                        interactions.add(new Edge<>(aminoAcid1, aminoAcid2, scheme.getWeight(plipInteraction))));
            }
        }

        return new Graph<>(aminoAcids, interactions);
    }

    public static PartitionedGraph<AminoAcid> createPartitionedGraphFromSecondaryStructureElements(Chain chain) {
        logger.info("partitioning {} using secondary structure elements",
                chain.getChainIdentifier());
        // annotate using DSSP
        new DictionaryOfProteinSecondaryStructure().process(chain.getParentStructure());

        List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());
        List<Module<AminoAcid>> modules = new ArrayList<>();
        List<AminoAcid> currentModule = new ArrayList<>();
        String lastSse = null;

        for(AminoAcid aminoAcid : aminoAcids) {
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
                    modules.add(new Module<>(String.valueOf(modules.size() + 1), currentModule));
                }
                currentModule = new ArrayList<>();
                currentModule.add(aminoAcid);
            }
        }

        modules.add(new Module<>(String.valueOf(modules.size() + 1), currentModule));

        // for each module: connect all contained nodes
        List<Edge<AminoAcid>> edges = createEdgesNaively(modules);

        return new PartitionedGraph<>(new Graph<>(aminoAcids, edges), modules);
    }

    private static List<Edge<AminoAcid>> createEdgesNaively(List<Module<AminoAcid>> modules) {
        return modules.stream()
                .flatMap(module -> SetOperations.uniquePairsOf(module.getNodes())
                        .map(pair -> new Edge<>(pair.getLeft(), pair.getRight())))
                .collect(Collectors.toList());
    }

    public static PartitionedGraph<AminoAcid> createPartitionedGraphFromPlipData(Chain chain,
                                                                                 String plipDocument,
                                                                                 double inflation,
                                                                                 WeightingScheme scheme) {
        logger.info("partitioning {} using MCL on PLIP interactions with scheme {}",
                chain.getChainIdentifier(),
                scheme);
        return createPartitionedGraphFromPlip(createGraphFromPlipDocument(chain, plipDocument, scheme), inflation);
    }

    public static PartitionedGraph<AminoAcid> createPartitionedGraphFromPlipData(Chain chain,
                                                                                 Path plipFile,
                                                                                 double inflation,
                                                                                 WeightingScheme scheme) {
        try {
            String plipDocument = Files.lines(plipFile).collect(Collectors.joining(System.lineSeparator()));
            return createPartitionedGraphFromPlipData(chain, plipDocument, inflation, scheme);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    private static PartitionedGraph<AminoAcid> createPartitionedGraphFromPlip(Graph<AminoAcid> graph, double inflation) {
        MCL mcl = new MCL(MCL.DEFAULT_EXPAND_FACTOR,
                inflation,
                MCL.DEFAULT_MULT_FACTOR,
                MCL.DEFAULT_MAX_ITERATIONS,
                MCL.DEFAULT_SIMILARITY_FUNCTION);

        return mcl.partitionGraph(graph);
    }


    public static PartitionedGraph<AminoAcid> createPartitionedGraphFromDefinitionFile(Chain chain, Path definitionFile) {
        logger.info("partitioning {} using definition file {}",
                chain.getChainIdentifier(),
                definitionFile);
        try {
            List<String> lines = Files.readAllLines(definitionFile);
            List<AminoAcid> aminoAcids = new ArrayList<>();
            List<Module<AminoAcid>> modules = new ArrayList<>();
            for(String range : lines) {
                if(range.startsWith("#")) {
                    continue;
                }

                String id = range.split(":")[0];
                String rawRanges = range.split(":")[1];

                List<AminoAcid> nodes = Pattern.compile(",").splitAsStream(rawRanges)
                        .map(rawRange -> rawRange.split("-"))
                        .flatMap(rawRange -> IntStream.range(Integer.valueOf(rawRange[0]), Integer.valueOf(rawRange[1]) + 1).boxed())
                        .map(residueNumber -> chain.select().residueNumber(residueNumber).asAminoAcid())
                        .collect(Collectors.toList());
                aminoAcids.addAll(nodes);
                modules.add(new Module<>(id, nodes));
            }

            List<Edge<AminoAcid>> edges = createEdgesNaively(modules);

            return new PartitionedGraph<>(new Graph<>(aminoAcids, edges), modules);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }
}
