package de.bioforscher.jstructure.mutation;

import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;

import java.util.stream.Collectors;

/**
 * The environment around a particular residue.
 * Created by bittrich on 7/6/17.
 */
public class Environment {
    private final AminoAcid centralGroup;
    private final GroupContainer sequentialNeighbors;
    private final String sequence;
    private final GroupContainer spatialNeighbors;
    private final PLIPInteractionContainer interactionNeighbors;

    public Environment(AminoAcid centralGroup,
                       GroupContainer sequentialNeighbors,
                       String sequence,
                       GroupContainer spatialNeighbors,
                       PLIPInteractionContainer interactionNeighbors) {
        this.centralGroup = centralGroup;
        this.sequentialNeighbors = sequentialNeighbors;
        this.sequence = sequence;
        this.spatialNeighbors = spatialNeighbors;
        this.interactionNeighbors = interactionNeighbors;
    }

    public AminoAcid getCentralGroup() {
        return centralGroup;
    }

    public GroupContainer getSequentialNeighbors() {
        return sequentialNeighbors;
    }

    public String getSequence() {
        return sequence;
    }

    public GroupContainer getSpatialNeighbors() {
        return spatialNeighbors;
    }

    public PLIPInteractionContainer getInteractionNeighbors() {
        return interactionNeighbors;
    }

    @Override
    public String toString() {
        return "Environment{" +
                "centralGroup=" + centralGroup +
                ",\n sequentialNeighbors=" + sequentialNeighbors.groups().map(Object::toString).collect(Collectors.joining(", ", "[", "]")) +
                ",\n sequence='" + sequence + '\'' +
                ",\n spatialNeighbors=" + spatialNeighbors.groups().map(Object::toString).collect(Collectors.joining(", ", "[", "]")) +
                ",\n interactionNeighbors=" + interactionNeighbors.getInteractions().stream().map(Object::toString).collect(Collectors.joining(", ", "[", "]")) +
                '}';
    }
}
