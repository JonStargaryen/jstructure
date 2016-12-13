package design.aggregation;

import de.bioforscher.jstructure.alignment.AlignmentResult;
import de.bioforscher.jstructure.alignment.svd.SVDSuperimposer;
import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.parser.ProteinParser;
import design.DesignConstants;

import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.List;
import java.util.Locale;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Aligns fragments against its cluster consensus fragment (if any).
 * Created by S on 09.12.2016.
 */
public class S09_AlignFragmentsByConsensus {
    public static void main(String[] args) {
        S06_ExtractSequences.TOPOLOGIES.stream()
                .limit(1)
                .forEach(S09_AlignFragmentsByConsensus::handleTopology);
    }

    private static final SVDSuperimposer SVD_SUPERIMPOSER = new SVDSuperimposer();
    private static final DecimalFormat DECIMAL_FORMAT = new DecimalFormat("0.000", DecimalFormatSymbols.getInstance(Locale.US));

    private static void handleTopology(final String topology) {
        System.out.println("topology: " + topology);
        Stream.of(SequenceMotifDefinition.values())
                .map(SequenceMotifDefinition::name)
                .forEach(motif -> {
                        List<AtomContainer> consensusFragments = loadConsensusFragments(topology, motif);
                        System.out.println("motif: " + motif);
                        System.out.println("aligning structures");

                        if(consensusFragments.size() == 0) {
                            return;
                        }

                        // align fragment to the best fitting consensus
                        List<AlignmentResult> alignedProteins = DesignConstants.list(Paths.get(DesignConstants.MOTIF_FRAGMENT_BY_TOPOLOGY_DIR +
                                topology + "/"))
                                .filter(path -> path.getFileName().toString().startsWith(motif))
                                .map(ProteinParser::parsePDBFile)
                                .filter(protein -> protein.getSize() != 0)
                                // align each to its most similar consensus fragment
                                .map(protein -> consensusFragments.stream()
                                        .map(consensus -> SVD_SUPERIMPOSER.align(consensus, protein))
                                        .reduce((a1, a2) -> a1.getRmsd() < a2.getRmsd() ? a1 : a2)
                                        .get())
                                .collect(Collectors.toList());

                        System.out.println("writing aligned files");
                        // write structures
                        alignedProteins.stream()
                                .filter(alignmentResult -> alignmentResult.getRmsd() < 1)
                                .forEach(alignmentResult -> DesignConstants.write(Paths.get(DesignConstants.ALIGNED_MOTIF_FRAGMENT_BY_TOPOLOGY_DIR +
                                        topology + "/" + alignmentResult.getOriginalQuery().getIdentifier() + "-" +
                                        alignmentResult.getOriginalReference().getIdentifier().split("-")[2] + "-" + DECIMAL_FORMAT.format(alignmentResult.getRmsd())
                                        + DesignConstants.PDB_SUFFIX),
                                        alignmentResult.getAlignedQuery().composePDBRecord().getBytes()));
                });
    }

    private static List<AtomContainer> loadConsensusFragments(String topology, String motif) {
        List<AtomContainer> fragments = DesignConstants.list(Paths.get(DesignConstants.ALIGNED_MOTIF_FRAGMENT_CONSENSUS_DIR + topology + "/"))
                .filter(path -> path.toFile().getName().startsWith(motif))
                .filter(path -> !path.toFile().getName().split("-")[1].equals("0"))
                .filter(path -> Integer.valueOf(path.toFile().getName().split("-")[2].split("\\.")[0]) >= DesignConstants.RARE_CLUSTER_THRESHOLD)
                .map(ProteinParser::parsePDBFile)
                .collect(Collectors.toList());

        if(fragments.size() < 2) {
            return fragments;
        }

        AtomContainer reference = fragments.get(0);
        return fragments.stream()
                .map(fragment -> SVD_SUPERIMPOSER.align(reference, fragment))
                .map(AlignmentResult::getAlignedQuery)
                .collect(Collectors.toList());
    }
}
