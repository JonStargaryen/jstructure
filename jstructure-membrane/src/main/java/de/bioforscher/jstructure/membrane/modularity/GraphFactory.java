package de.bioforscher.jstructure.membrane.modularity;

import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.interactions.PLIPIntraMolecularAnnotator;
import de.bioforscher.jstructure.feature.sse.SecondaryStructureElement;
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
import java.util.Collections;
import java.util.List;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Constructs graph representations of protein chains.
 */
public class GraphFactory {
    private static final Logger logger = LoggerFactory.getLogger(GraphFactory.class);
    private static final PLIPIntraMolecularAnnotator plipIntraMolecularAnnotator = new PLIPIntraMolecularAnnotator();

    public static Graph<AminoAcid> createGraphFromPlipDocument(Chain chain, String plipDocument) {
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

                if(!Collections.disjoint(interactions1.getInteractions(), interactions2.getInteractions())) {
                    interactions.add(new Edge<>(aminoAcid1, aminoAcid2));
                }
            }
        }

        return new Graph<>(aminoAcids, interactions);
    }

    public static PartitionedGraph<AminoAcid> createPartitionedGraphFromSecondaryStructureElements(Chain chain) {
        // annotate using DSSP
        new DictionaryOfProteinSecondaryStructure().process(chain.getParentStructure());

        List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());
        List<Module<AminoAcid>> modules = new ArrayList<>();
        List<AminoAcid> currentModule = new ArrayList<>();
        SecondaryStructureElement lastSse = null;

        for(AminoAcid aminoAcid : aminoAcids) {
            SecondaryStructureElement sse = aminoAcid.getFeature(DSSPSecondaryStructure.class).getSecondaryStructure();

            // we define modules as change of sse
            if(sse.equals(lastSse)) {
                // when it matches, append last module
                currentModule.add(aminoAcid);
            } else {
                // if not: change sse, assign module and clear list
                lastSse = sse;
                modules.add(new Module<>(String.valueOf(modules.size() + 1), new Graph<>(currentModule, new ArrayList<>())));
                currentModule = new ArrayList<>();
            }
        }

        return new PartitionedGraph<>(new Graph<>(aminoAcids, new ArrayList<>()), modules);
    }

    public static PartitionedGraph<AminoAcid> createPartitionedGraphFromMCL(Chain chain, Path plipFile) {
        return createPartitionedGraphFromMCL(chain, plipFile, MCL.DEFAULT_INFLATE_FACTOR);
    }

    public static PartitionedGraph<AminoAcid> createPartitionedGraphFromMCL(Chain chain, Path plipFile, double inflation) {
        try {
            String plipDocument = Files.lines(plipFile).collect(Collectors.joining(System.lineSeparator()));
            Graph<AminoAcid> graph = createGraphFromPlipDocument(chain, plipDocument);

            MCL mcl = new MCL(MCL.DEFAULT_EXPAND_FACTOR,
                    inflation,
                    MCL.DEFAULT_MULT_FACTOR,
                    MCL.DEFAULT_MAX_ITERATIONS,
                    MCL.DEFAULT_SIMILARITY_FUNCTION);

            return mcl.clusterGraph(graph);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

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

            return new PartitionedGraph<>(new Graph<>(aminoAcids, new ArrayList<>()), modules);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    public static PartitionedGraph<AminoAcid> createPartitionedGraphFromDefinitionFile(Chain chain, Path definitionFile) {
        try {
            List<String> lines = Files.readAllLines(definitionFile);
            List<AminoAcid> aminoAcids = new ArrayList<>();
            List<Module<AminoAcid>> modules = new ArrayList<>();
            for(String range : lines) {
                String id = range.split(":")[0];
                String rawRanges = range.split(":")[1];

                List<AminoAcid> nodes = Pattern.compile(",").splitAsStream(rawRanges)
                        .map(rawRange -> rawRange.split("-"))
                        .flatMap(rawRange -> IntStream.range(Integer.valueOf(rawRange[0]), Integer.valueOf(rawRange[1]) + 1).boxed())
                        .map(residueNumber -> chain.select().residueNumber(residueNumber).asAminoAcid())
                        .collect(Collectors.toList());
                aminoAcids.addAll(nodes);
                modules.add(new Module<>(id, new Graph<>(nodes, new ArrayList<>())));
            }

            return new PartitionedGraph<>(new Graph<>(aminoAcids, new ArrayList<>()), modules);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }
}
