package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.structure.filter.ResiduePairDistanceCutoffFilter;
import de.bioforscher.jstructure.model.structure.scheme.AlphaCarbonRepresentationScheme;
import de.bioforscher.jstructure.model.structure.scheme.RepresentationScheme;

import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Specifies the capabilities of a group container (mostly a {@link Chain}).
 * Created by S on 30.09.2016.
 */
public interface GroupContainer extends AtomContainer {
    /**
     * Access to all groups associated to this container.
     * @return a stream of groups
     */
    Stream<Group> groups();

    Stream<Residue> residues();

    /**
     *
     * @param residueNumber
     * @return
     */
    default Optional<Residue> residue(int residueNumber) {
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
     * Returns a stream of all residue combination whose euclidean distance is below a defined threshold.
     * @param distanceCutoff the distance threshold below two atoms are considered to be in contact
     * @param representationScheme how to represent groups?
     * @return all residue pairs in contact according to the distance cutoff
     * @see Pair#uniquePairsOf(List)
     * @see de.bioforscher.jstructure.model.structure.filter.GroupPairDistanceCutoffFilter
     */
    default Stream<Pair<Residue, Residue>> residuePairsInContact(double distanceCutoff,
                                                                 RepresentationScheme representationScheme) {
        final ResiduePairDistanceCutoffFilter distanceCutoffFilter = new ResiduePairDistanceCutoffFilter(distanceCutoff,
                representationScheme);
        return Pair.uniquePairsOf(residues().collect(Collectors.toList())).filter(distanceCutoffFilter);
    }

    /**
     * @see GroupContainer#residuePairsInContact(double, RepresentationScheme)
     */
    default Stream<Pair<Residue, Residue>> residuePairsInContact(double distanceCutoff) {
        return residuePairsInContact(distanceCutoff, new AlphaCarbonRepresentationScheme());
    }

    /**
     * Returns all residue pairs.
     * @return a stream of all residue pairs
     */
    default Stream<Pair<Residue, Residue>> residuePairs() {
        return Pair.uniquePairsOf(residues().collect(Collectors.toList()));
    }
}