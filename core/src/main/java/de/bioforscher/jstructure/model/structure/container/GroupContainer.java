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
        return groups()
                .filter(Group::isAminoAcid)
                .map(Group::getGroupInformation)
                .map(GroupInformation::getOneLetterCode)
                .collect(Collectors.joining());
    }

    @Override
    default GroupContainer getCopy() {
        return (GroupContainer) AtomContainer.super.getCopy();
    }
}