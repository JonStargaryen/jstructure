package de.bioforscher.jstructure.model.structure.container;

import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.structure.AminoAcid;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Residue;
import de.bioforscher.jstructure.model.structure.filter.ResiduePairDistanceCutoffFilter;
import de.bioforscher.jstructure.model.structure.scheme.AlphaCarbonRepresentationScheme;
import de.bioforscher.jstructure.model.structure.scheme.RepresentationScheme;

import java.util.List;
import java.util.NoSuchElementException;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Specifies the capabilities of a group container (mostly a {@link Chain}).
 * Created by S on 30.09.2016.
 */
public interface GroupContainer extends AtomContainer {
    List<Group> getGroups();

    default Stream<Residue> residues() {
        return groups().filter(Residue.class::isInstance)
                       .map(Residue.class::cast)
                       //TODO replace this once the parser is working
                       .filter(residue -> residue.getAminoAcid().isStandardAminoAcid());
    }

    default List<Residue> getResidues() {
        return residues().collect(Collectors.toList());
    }

    /**
     * Access to all groups associated to this container.
     * @return a stream of groups
     */
    default Stream<Group> groups() {
        return getGroups().stream();
    }

    /**
     *
     * @param residueNumber
     * @return
     */
    default Residue getResidue(final int residueNumber) {
        return findResidue(residueNumber).orElseThrow(NoSuchElementException::new);
    }

    default Optional<Residue> findResidue(final int residueNumber) {
        return residues().filter(residue -> residue.getResidueNumber() == residueNumber).findFirst();
    }

    /**
     *
     * @return
     */
    default String getSequence() {
        return residues().map(Residue::getAminoAcid)
                         .map(AminoAcid::getOneLetterCode)
                         .collect(Collectors.joining());
    }

    /**
     * Returns a stream of all getResidue combination whose euclidean distance is below a defined threshold.
     * @param distanceCutoff the distance threshold below two atoms are considered to be in contact
     * @param representationScheme how to represent groups?
     * @return all getResidue pairs in contact according to the distance cutoff
     * @see Pair#uniquePairsOf(List)
     * @see de.bioforscher.jstructure.model.structure.filter.GroupPairDistanceCutoffFilter
     */
    default Stream<Pair<Residue, Residue>> residuePairsInContact(double distanceCutoff,
                                                                 RepresentationScheme representationScheme) {
        final ResiduePairDistanceCutoffFilter distanceCutoffFilter =
                new ResiduePairDistanceCutoffFilter(distanceCutoff, representationScheme);
        return residuePairs().filter(distanceCutoffFilter);
    }

    /**
     * @see GroupContainer#residuePairsInContact(double, RepresentationScheme)
     */
    default Stream<Pair<Residue, Residue>> residuePairsInContact(double distanceCutoff) {
        return residuePairsInContact(distanceCutoff, new AlphaCarbonRepresentationScheme());
    }

    /**
     * Returns all getResidue pairs.
     * @return a stream of all getResidue pairs
     */
    default Stream<Pair<Residue, Residue>> residuePairs() {
        return Pair.uniquePairsOf(getResidues());
    }
}