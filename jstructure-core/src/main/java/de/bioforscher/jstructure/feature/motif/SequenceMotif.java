package de.bioforscher.jstructure.feature.motif;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.StructureCollectors;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.family.GroupInformation;

import java.util.stream.Collectors;

/**
 * The container object of a found sequence motif.
 * Created by S on 02.10.2016.
 */
public class SequenceMotif {
    private final SequenceMotifDefinition motifDefinition;
    private final Group startGroup;
    private final Group endGroup;
    private final GroupContainer groups;

    public SequenceMotif(SequenceMotifDefinition candidate, Group startGroup, Group endGroup) {
        this.motifDefinition = candidate;
        this.startGroup = startGroup;
        this.endGroup = endGroup;

        int startResNum = startGroup.getResidueNumber();
        int endResNum = endGroup.getResidueNumber();
        Chain chain = startGroup.getParentChain();

        // extract groups
        this.groups = chain.groups()
                .filter(group -> group.getResidueNumber() >= startResNum && group.getResidueNumber() <= endResNum)
                .collect(StructureCollectors.toGroupContainer());
    }

    public GroupContainer getGroupContainer() {
        return groups;
    }

    public String getSequence() {
        return groups.groups()
                .map(Group::getGroupInformation)
                .map(GroupInformation::getOneLetterCode)
                .collect(Collectors.joining());
    }

    public Group getEndGroup() {
        return endGroup;
    }

    public Group getStartGroup() {
        return startGroup;
    }

    public SequenceMotifDefinition getMotifDefinition() {
        return motifDefinition;
    }

    @Override
    public String toString() {
        return getClass().getSimpleName() + " motifDefinition='" + motifDefinition + "' sequence='" + getSequence() +
                "' startGroup='" + startGroup.getThreeLetterCode() + "-" + startGroup.getResidueNumber() +
                "' endGroup='" + endGroup.getThreeLetterCode() + "-" + endGroup.getResidueNumber() + "'";
    }
}