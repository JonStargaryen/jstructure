package de.bioforscher.jstructure.membrane.modularity;

import de.bioforscher.jstructure.feature.interactions.PLIPInteraction;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.interactions.PLIPIntraMolecularAnnotator;
import de.bioforscher.jstructure.feature.sse.dssp.DSSPSecondaryStructure;
import de.bioforscher.jstructure.feature.sse.dssp.DictionaryOfProteinSecondaryStructure;
import de.bioforscher.jstructure.mathematics.graph.Edge;
import de.bioforscher.jstructure.mathematics.graph.Graph;
import de.bioforscher.jstructure.mathematics.graph.PartitionedGraph;
import de.bioforscher.jstructure.mathematics.graph.clustering.Module;
import de.bioforscher.jstructure.mathematics.graph.clustering.algorithms.MCL;
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
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Constructs graph representations of protein chains.
 */
public class GraphFactory {
    //TODO value
    private static final double COVALENT_BOND_ENERGY_CONTRIBUTION = 50;
    private static final Logger logger = LoggerFactory.getLogger(GraphFactory.class);
    private static final PLIPIntraMolecularAnnotator plipIntraMolecularAnnotator = new PLIPIntraMolecularAnnotator();

    public static Graph<AminoAcid> createUnweightedGraphFromPlipDocument(Chain chain, String plipDocument) {
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
                    interactions.add(new Edge<>(aminoAcid1, aminoAcid2));
                    continue;
                }

                PLIPInteractionContainer interactions2 = aminoAcid2.getFeature(PLIPInteractionContainer.class);

                Optional<PLIPInteraction> potentialInteraction = interactions1.getInteractions()
                        .stream()
                        .filter(plipInteraction -> interactions2.getInteractions().contains(plipInteraction))
                        .findFirst();
                potentialInteraction.ifPresent(plipInteraction ->
                        interactions.add(new Edge<>(aminoAcid1, aminoAcid2)));
            }
        }

        return new Graph<>(aminoAcids, interactions);
    }

    public static Graph<AminoAcid> createWeightedGraphFromPlipDocument(Chain chain, String plipDocument) {
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
                    interactions.add(new Edge<>(aminoAcid1, aminoAcid2, COVALENT_BOND_ENERGY_CONTRIBUTION));
                    continue;
                }

                PLIPInteractionContainer interactions2 = aminoAcid2.getFeature(PLIPInteractionContainer.class);

                Optional<PLIPInteraction> potentialInteraction = interactions1.getInteractions()
                        .stream()
                        .filter(plipInteraction -> interactions2.getInteractions().contains(plipInteraction))
                        .findFirst();
                potentialInteraction.ifPresent(plipInteraction ->
                        interactions.add(new Edge<>(aminoAcid1, aminoAcid2, plipInteraction.getEnergyContribution())));
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

        return new PartitionedGraph<>(aminoAcids, modules);
    }

    public static PartitionedGraph<AminoAcid> createPartitionedGraphFromUnweightedPlipData(Chain chain, Path plipFile, double inflation) {
        logger.info("partitioning {} using MCL on weighted PLIP interactions",
                chain.getChainIdentifier());
        try {
            String plipDocument = Files.lines(plipFile).collect(Collectors.joining(System.lineSeparator()));
            return createPartitionedGraphFromPlip(createUnweightedGraphFromPlipDocument(chain, plipDocument), inflation);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    public static PartitionedGraph<AminoAcid> createPartionedGraphFromWeightedPlipData(Chain chain, Path plipFile, double inflation) {
        logger.info("partitioning {} using MCL on weighted PLIP interactions",
                chain.getChainIdentifier());
        try {
            String plipDocument = Files.lines(plipFile).collect(Collectors.joining(System.lineSeparator()));
            return createPartitionedGraphFromPlip(createWeightedGraphFromPlipDocument(chain, plipDocument), inflation);
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

        return mcl.clusterGraph(graph);
    }

    @Deprecated
    public static PartitionedGraph<AminoAcid> createPartitionedGraphFromNetCartoFile(Chain chain, Path netCartoFile) {
        try {
            List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());
            List<Module<AminoAcid>> modules = Files.lines(netCartoFile)
                    .filter(line -> !line.startsWith("#"))
                    .map(line -> new Module<>(line.split(" ")[0],
                            new Graph<>(Pattern.compile("\\s+").splitAsStream(line.split("---")[1].trim())
                                    .mapToInt(Integer::valueOf)
                                    .mapToObj(residueNumber -> chain.select()
                                            .residueNumber(residueNumber)
                                            .asAminoAcid())
                                    .collect(Collectors.toList()), new ArrayList<>())))
                    .collect(Collectors.toList());

            return new PartitionedGraph<>(aminoAcids, modules);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
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

            return new PartitionedGraph<>(aminoAcids, modules);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }
}
