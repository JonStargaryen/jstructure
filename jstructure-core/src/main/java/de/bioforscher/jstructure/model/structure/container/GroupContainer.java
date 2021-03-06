package de.bioforscher.jstructure.model.structure.container;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Water;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.nucleotide.Nucleotide;

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
    default Stream<AminoAcid> aminoAcids() {
        return groups()
                .filter(Group::isAminoAcid)
                .map(AminoAcid.class::cast);
    }

    default Stream<Group> ligands() {
        return groups()
                .filter(Group::isLigand);
    }

    default Stream<Nucleotide> nucleotides() {
        return groups()
                .filter(Group::isNucleotide)
                .map(Nucleotide.class::cast);
    }

    default Stream<Water> waters() {
        return groups()
                .filter(Group::isWater)
                .map(Water.class::cast);
    }

    default String getAminoAcidSequence() {
        return aminoAcids()
                .map(AminoAcid::getOneLetterCode)
                .collect(Collectors.joining());
    }

    default String getNucleotideSequence() {
        return nucleotides()
                .map(Nucleotide::getOneLetterCode)
                .collect(Collectors.joining());
    }
}