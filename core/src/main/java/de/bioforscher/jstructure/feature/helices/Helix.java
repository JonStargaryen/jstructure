package de.bioforscher.jstructure.feature.helices;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.StructureCollectors;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.family.GroupInformation;

import java.util.stream.Collectors;

/**
 * Properties of a helix.
 * Created by bittrich on 2/13/17.
 */
public class Helix {
    private final Group startGroup;
    private final Group endGroup;
    private final GroupContainer groupContainer;

    public Helix(Group startGroup, Group endGroup) {
        this.startGroup = startGroup;
        this.endGroup = endGroup;

        int startResNum = startGroup.getResidueNumber();
        int endResNum = endGroup.getResidueNumber();
        Chain chain = startGroup.getParentChain();

        // extract groups
        this.groupContainer = chain.groups()
                .filter(group -> group.getResidueNumber() >= startResNum && group.getResidueNumber() <= endResNum)
                .collect(StructureCollectors.toGroupContainer());
    }

    public GroupContainer getGroupContainer() {
        return groupContainer;
    }

    public String getSequence() {
        return groupContainer.groups()
                .map(Group::getGroupInformation)
                .map(GroupInformation::getOneLetterCode)
                .collect(Collectors.joining());
    }

    public Group getStartGroup() {
        return startGroup;
    }

    public Group getEndGroup() {
        return endGroup;
    }

    @Override
    public String toString() {
        return getClass().getSimpleName() + " sequence='" + getSequence() + "' startGroup='" +
                startGroup.getIdentifier() + "' endGroup='" + endGroup.getIdentifier() + "'";
    }
}
