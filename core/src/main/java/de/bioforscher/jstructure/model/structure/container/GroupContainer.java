package de.bioforscher.jstructure.model.structure.container;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.family.GroupInformation;
import de.bioforscher.jstructure.model.structure.selection.Selection;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Specifies the capabilities of a group container (mostly a {@link Chain}).
 * Created by S on 30.09.2016.
 */
public interface GroupContainer extends AtomContainer {
    /**
     * Never manipulate the returned collection as it is not guaranteed the actually modify the internal list(s).
     * @return all associated groups
     */
    List<Group> getGroups();

    /**
     * Access to all groups associated to this container.
     * @return all groups
     */
    default Stream<Group> groups() {
        return getGroups().stream();
    }

    /**
     * Access to all amino acids associated to this container
     * @return all amino acids
     */
    default Stream<Group> aminoAcids() {
        return Selection.on(this).aminoAcids().asFilteredGroups();
    }

    default String getAminoAcidSequence() {
        return aminoAcids()
                .map(Group::getGroupInformation)
                .map(GroupInformation::getOneLetterCode)
                .collect(Collectors.joining());
    }

    @Override
    default GroupContainer getCopy() {
        return (GroupContainer) AtomContainer.super.getCopy();
    }
}