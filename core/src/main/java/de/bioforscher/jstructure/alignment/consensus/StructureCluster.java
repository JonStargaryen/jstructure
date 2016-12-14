package de.bioforscher.jstructure.alignment.consensus;

import de.bioforscher.jstructure.alignment.AlignmentResult;
import de.bioforscher.jstructure.alignment.svd.SVDSuperimposer;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

/**
 * Represents a structure cluster created by {@link FragmentClusteringComposer}.
 * Created by S on 05.12.2016.
 */
public class StructureCluster {
    private List<AtomContainer> entries;
    private AtomContainer consensus;
    private SVDSuperimposer svdSuperimposer;

    public static final StructureCluster EMPTY_CLUSTER = new StructureCluster();

    private StructureCluster() {
        this.entries = new ArrayList<>();
        this.svdSuperimposer = new SVDSuperimposer(true);
    }

    StructureCluster(AtomContainer container) {
        this();
        this.entries.add(Objects.requireNonNull(container));
        this.consensus = container;
    }

    public AtomContainer add(AtomContainer container) {
        entries.add(container);
        consensus = updateConsensus(container);
        return consensus;
    }

    public AtomContainer getConsensusRepresentation() {
        return consensus;
    }

    public List<AtomContainer> getEntries() {
        return entries;
    }

    private AtomContainer updateConsensus(AtomContainer container) {
        //TODO we should store fragments aligned relatively to the consensus, so it is not aligned twice
        AlignmentResult alignment = svdSuperimposer.align(consensus, container);
        // merge new entry to existing consensus
        //TODO this may be fast, but potentially merging everything is more robust
        //TODO delegating to the abstract class is not really nice
        return AbstractConsensusComposer.mergeContainerPair(alignment.getOriginalReference(), alignment.getOriginalQuery());
    }

    double getMinimalDistance(AtomContainer container) {
        return entries.stream()
                .map(entry -> svdSuperimposer.align(entry, container))
                .mapToDouble(AlignmentResult::getRmsd)
                .min()
                .orElseThrow(() -> new IllegalArgumentException("rmsd to clusters could not be computed"));
    }

    double getDistanceToConsensus(AtomContainer container) {
        return svdSuperimposer.align(consensus, container).getRmsd();
    }
}