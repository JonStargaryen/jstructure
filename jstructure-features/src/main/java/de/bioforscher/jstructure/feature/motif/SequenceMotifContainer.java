package de.bioforscher.jstructure.feature.motif;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.identifier.ChainIdentifier;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * A collection of all {@link SequenceMotif} instances associated to a given instance.
 * Created by bittrich on 5/17/17.
 */
public class SequenceMotifContainer extends FeatureContainerEntry {
    private List<SequenceMotif> sequenceMotifs;

    private SequenceMotifContainer(AbstractFeatureProvider featureProvider, List<SequenceMotif> sequenceMotifs) {
        super(featureProvider);
        this.sequenceMotifs = sequenceMotifs;
    }

    SequenceMotifContainer(AbstractFeatureProvider featureProvider) {
        super(featureProvider);
        this.sequenceMotifs = new ArrayList<>();
    }

    public List<SequenceMotif> getSequenceMotifs() {
        return sequenceMotifs;
    }

    void addSequenceMotif(SequenceMotif sequenceMotif) {
        sequenceMotifs.add(sequenceMotif);
    }

    boolean containsSequenceMotif(SequenceMotif sequenceMotif) {
        return sequenceMotifs.contains(sequenceMotif);
    }

    public SequenceMotifContainer getEmbeddingSequenceMotifsFor(Group group) {
        ChainIdentifier chainId = group.getParentChain().getChainIdentifier();
        int residueNumber = group.getResidueIdentifier().getResidueNumber();
        List<SequenceMotif> filteredSequenceMotifs = sequenceMotifs.stream()
                // filter for correct chain
                .filter(sequenceMotif -> sequenceMotif.getChainId().equals(chainId))
                // filter for correct residue number range
                .filter(sequenceMotif -> sequenceMotif.getStartResidueNumber() <= residueNumber &&
                        sequenceMotif.getEndResidueNumber() >= residueNumber)
                .collect(Collectors.toList());
        return new SequenceMotifContainer(getFeatureProvider(), filteredSequenceMotifs);
    }
}
