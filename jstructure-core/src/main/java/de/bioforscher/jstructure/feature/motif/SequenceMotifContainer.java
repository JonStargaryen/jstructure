package de.bioforscher.jstructure.feature.motif;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;

import java.util.ArrayList;
import java.util.List;

/**
 * A collection of all {@link SequenceMotif} instances associated to a given instance.
 * Created by bittrich on 5/17/17.
 */
public class SequenceMotifContainer extends FeatureContainerEntry {
    private List<SequenceMotif> sequenceMotifs;

    public SequenceMotifContainer(AbstractFeatureProvider featureProvider) {
        super(featureProvider);
        this.sequenceMotifs = new ArrayList<>();
    }

    public List<SequenceMotif> getSequenceMotifs() {
        return sequenceMotifs;
    }

    public void addSequenceMotif(SequenceMotif sequenceMotif) {
        sequenceMotifs.add(sequenceMotif);
    }

    public boolean containsSequenceMotif(SequenceMotif sequenceMotif) {
        return sequenceMotifs.contains(sequenceMotif);
    }

    //TODO convenience functions to check if some group lies in any (or specific) motif etc
}
