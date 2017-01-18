package de.bioforscher.jstructure.model.structure.selection;

import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.model.Combinatorics;
import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.structure.*;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.container.ChainContainer;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.family.AminoAcidFamily;
import de.bioforscher.jstructure.model.structure.scheme.AlphaCarbonRepresentationScheme;
import de.bioforscher.jstructure.model.structure.scheme.BetaCarbonRepresentationScheme;
import de.bioforscher.jstructure.model.structure.scheme.RepresentationScheme;

import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * The specification and implementation of the selection API.
 * Created by S on 21.11.2016.
 */
public class Selection {
    public static AtomSelection on(AtomContainer atomContainer) {
        return new AtomSelection(Objects.requireNonNull(atomContainer));
    }

    public static AtomPairSelection pairsOn(AtomContainer atomContainer) {
        return new AtomPairSelection(Objects.requireNonNull(atomContainer));
    }

    public static GroupSelection on(GroupContainer groupContainer) {
        return new GroupSelection(Objects.requireNonNull(groupContainer));
    }

    public static GroupPairSelection pairsOn(GroupContainer groupContainer) {
        return new GroupPairSelection(Objects.requireNonNull(groupContainer));
    }

    public static ChainSelection on(ChainContainer chainContainer) {
        return new ChainSelection(Objects.requireNonNull(chainContainer));
    }

    /**
     * Merges a collection of predicates by chaining them via {@link Predicate#and(Predicate)}.
     * @param predicates a number of predicates
     * @param <T> the type these predicates operate upon
     * @return a predicate which needs all predicates to evaluate to <code>true</code> - if no predicates were provided,
     *      the expression will evaluate to <code>true</code> by default
     */
    private static <T> Predicate<T> compose(List<Predicate<T>> predicates) {
        return predicates.stream()
                .reduce(Predicate::and)
                .orElse(x -> true);
    }

    public static class AtomSelection {
        AtomContainer atomContainer;
        GroupContainer groupContainer;
        ChainContainer chainContainer;
        List<Predicate<Atom>> atomPredicates;
        List<Predicate<Group>> groupPredicates;
        List<Predicate<Chain>> chainPredicates;
        boolean cloneRequested;

        AtomSelection(AtomContainer atomContainer) {
            this.atomContainer = atomContainer;
            this.atomPredicates = new ArrayList<>();
            this.groupPredicates = new ArrayList<>();
            this.chainPredicates = new ArrayList<>();
        }

        public AtomSelection cloneElements() {
            cloneRequested = true;
            return this;
        }

        public Stream<Atom> asFilteredAtoms() {
            //TODO resolve this back-and-forth pain
            Stream<Group> prefilteredGroupStream = null;
            Stream<Atom> prefilteredAtomStream = null;
            // pre-filter groups, when there are requirements on the chains
            if(this instanceof ChainSelection) {
                prefilteredGroupStream = chainContainer.chains().filter(Selection.compose(chainPredicates)).flatMap(Chain::groups);
            }

            if(this instanceof GroupSelection) {
                // group selection, but no chain selection
                if(prefilteredGroupStream == null) {
                    prefilteredGroupStream = groupContainer.groups();
                }

                prefilteredAtomStream = prefilteredGroupStream.filter(Selection.compose(groupPredicates)).flatMap(Group::atoms);
            }

            // fallback, when no higher container is present
            if(prefilteredAtomStream == null) {
                prefilteredAtomStream = atomContainer.atoms();
            }

            Stream<Atom> filteredAtomStream = prefilteredAtomStream
                    .filter(Selection.compose(atomPredicates));
            if(cloneRequested) {
                return filteredAtomStream
                        .map(Atom::new);
            } else {
                return filteredAtomStream;
            }
        }

        public AtomContainer asAtomContainer() {
            return asFilteredAtoms()
                    .collect(StructureCollectors.toAtomContainer());
        }

        public AtomPairSelection asAtomPairSelection() {
            return pairsOn(asAtomContainer());
        }

        public Optional<Atom> asOptionalAtom() {
            return asFilteredAtoms()
                    .findAny();
        }

        public Atom asAtom() {
            return asOptionalAtom()
                    .orElseThrow(NoSuchElementException::new);
        }

        public AtomSelection element(Collection<Element> elements) {
            return element(elements.toArray(new Element[elements.size()]));
        }

        public AtomSelection element(Element... elements) {
            atomPredicates.add(atom -> Stream.of(elements)
                    .anyMatch(element -> atom.getElement().equals(element)));
            return this;
        }

        public AtomSelection atomName(Collection<String> atomNames){
            return atomName(atomNames.toArray(new String[atomNames.size()]));
        }

        public AtomSelection atomName(String... atomNames) {
            atomPredicates.add(atom -> Stream.of(atomNames)
                    .anyMatch(atomName -> atomName.equals(atom.getName())));
            return this;
        }

        public AtomSelection backboneAtoms() {
            atomPredicates.add(atom -> AminoAcidFamily.ATOM_NAMES.BACKBONE_ATOM_NAMES.contains(atom.getName()));
            return this;
        }

        public AtomSelection sideChainAtoms() {
            atomPredicates.add(atom -> AminoAcidFamily.ATOM_NAMES.SIDECHAIN_ATOM_NAMES.contains(atom.getName()));
            return this;
        }

        public AtomSelection allAtoms() {
            atomPredicates.add(atom -> Stream.concat(AminoAcidFamily.ATOM_NAMES.BACKBONE_ATOM_NAMES.stream(),
                    AminoAcidFamily.ATOM_NAMES.SIDECHAIN_ATOM_NAMES.stream()).collect(Collectors.toList()).contains(atom.getName()));
            return this;
        }

        public AtomSelection alphaCarbonAtoms() {
            atomPredicates.add(atom -> AminoAcidFamily.ATOM_NAMES.CA_ATOM_NAME.equals(atom.getName()));
            return this;
        }

        public AtomSelection backboneCarbonAtoms() {
            atomPredicates.add(atom -> AminoAcidFamily.ATOM_NAMES.C_ATOM_NAME.equals(atom.getName()));
            return this;
        }

        public AtomSelection backboneNitrogenAtoms() {
            atomPredicates.add(atom -> AminoAcidFamily.ATOM_NAMES.N_ATOM_NAME.equals(atom.getName()));
            return this;
        }

        public AtomSelection backboneOxygenAtoms() {
            atomPredicates.add(atom -> AminoAcidFamily.ATOM_NAMES.O_ATOM_NAME.equals(atom.getName()));
            return this;
        }

        public AtomSelection betaCarbonAtoms() {
            atomPredicates.add(atom -> AminoAcidFamily.ATOM_NAMES.CB_ATOM_NAME.equals(atom.getName()));
            return this;
        }

        public AtomSelection hydrogenAtoms() {
            atomPredicates.add(atom -> AminoAcidFamily.ATOM_NAMES.H_ATOM_NAMES.contains(atom.getName()));
            return this;
        }

        public AtomSelection backboneHydrogen() {
//            backboneAtoms(); - 18/1/17 fix: will never be true as backbone atoms no longer container hydrogen names
            hydrogenAtoms();
            return this;
        }

        public AtomSelection nonHydrogenAtoms() {
            atomPredicates.add(atom -> !AminoAcidFamily.ATOM_NAMES.H_ATOM_NAMES.contains(atom.getName()));
            return this;
        }

        public AtomSelection pdbSerial(int... pdbSerials) {
            atomPredicates.add(atom -> IntStream.of(pdbSerials).boxed().collect(Collectors.toList())
                    .contains(atom.getPdbSerial()));
            return this;
        }

        public AtomSelection pdbSerial(Range... pdbSerialRanges) {
            atomPredicates.add(atom -> Stream.of(pdbSerialRanges)
                    .anyMatch(range -> range.getLeft() >= atom.getPdbSerial() && range.getRight() <= atom.getPdbSerial()));
            return this;
        }

        public AtomSelection distance(Atom referenceAtom, double cutoff) {
            atomPredicates.add(atom -> LinearAlgebra3D.distanceFast(referenceAtom.getCoordinates(),
                    atom.getCoordinates()) < cutoff * cutoff);
            return this;
        }
    }

    public static class GroupSelection extends AtomSelection {
        GroupSelection(GroupContainer groupContainer) {
            super(groupContainer);
            this.groupContainer = groupContainer;
        }

        public GroupSelection cloneElements() {
            cloneRequested = true;
            return this;
        }

        public Stream<Group> asFilteredGroups() {
            Stream<Group> prefilteredGroupStream;
            // pre-filter groups, when there are requirements on the chains
            if(!(this instanceof ChainSelection)) {
                prefilteredGroupStream = groupContainer.groups();
            } else {
                prefilteredGroupStream = chainContainer.chains().filter(Selection.compose(chainPredicates)).flatMap(Chain::groups);
            }

            Stream<Group> filteredGroupStream =  prefilteredGroupStream
                    .filter(Selection.compose(groupPredicates));
            if(cloneRequested) {
                return filteredGroupStream
                        .map(Group::new);
            } else {
                return filteredGroupStream;
            }
        }

        public GroupContainer asGroupContainer() {
            return asFilteredGroups()
                    .collect(StructureCollectors.toGroupContainer());
        }

        public GroupPairSelection asGroupPairSelection() {
            return pairsOn(asGroupContainer());
        }

        public Optional<Group> asOptionalGroup() {
            return asFilteredGroups()
                    .findAny();
        }

        public Group asGroup() {
            return asOptionalGroup()
                    .orElseThrow(NoSuchElementException::new);
        }

        public GroupSelection aminoAcids() {
            groupPredicates.add(Group::isAminoAcid);
            return this;
        }

        public GroupSelection aminoAcids(AminoAcidFamily... aminoAcids) {
            // pre-filter for group type amino acid
            aminoAcids();
            groupPredicates.add(group -> Stream.of(aminoAcids)
                    .map(AminoAcidFamily::getThreeLetterCode)
                    .anyMatch(threeLetterCode -> threeLetterCode.equals(group.getThreeLetterCode())));
            return this;
        }

        public GroupSelection nucleotides() {
            groupPredicates.add(Group::isNucleotide);
            return this;
        }

        public GroupSelection hetatms() {
            groupPredicates.add(Group::isLigand);
            return this;
        }

        public GroupSelection groupName(String... groupNames) {
            groupPredicates.add(group -> Stream.of(groupNames)
                    .anyMatch(groupName -> groupName.equals(group.getThreeLetterCode())));
            return this;
        }

        public GroupSelection residueNumber(int... residueNumbers) {
            groupPredicates.add(group -> IntStream.of(residueNumbers).boxed().collect(Collectors.toList())
                    .contains(group.getResidueNumber()));
            return this;
        }

        public GroupSelection residueNumber(Range... residueNumberRanges) {
            groupPredicates.add(group -> Stream.of(residueNumberRanges)
                    .anyMatch(range -> range.getLeft() >= group.getResidueNumber() && range.getRight() <= group.getResidueNumber()));
            return this;
        }
    }

    public static class ChainSelection extends GroupSelection {
        ChainSelection(ChainContainer chainContainer) {
            super(chainContainer);
            this.chainContainer = chainContainer;
        }

        public ChainSelection cloneElements() {
            cloneRequested = true;
            return this;
        }

        public Stream<Chain> asFilteredChains() {
            Stream<Chain> filteredChainStream = chainContainer.chains()
                    .filter(Selection.compose(chainPredicates));
            if(cloneRequested) {
                return filteredChainStream
                        .map(Chain::new);
            } else {
                return filteredChainStream;
            }
        }

        public ChainContainer asChainContainer() {
            return asFilteredChains()
                    .collect(StructureCollectors.toChainContainer());
        }

        public Optional<Chain> asOptionalChain() {
            return asFilteredChains()
                    .findAny();
        }

        public Chain asChain() {
            return asOptionalChain()
                    .orElseThrow(NoSuchElementException::new);
        }

        public ChainSelection chainName(String... chainNames) {
            chainPredicates.add(chain -> Stream.of(chainNames)
                    .anyMatch(chainName -> chainName.equals(chain.getChainId())));
            return this;
        }
    }

    public static class AtomPairSelection {
        private List<Pair<Atom, Atom>> atomPairs;
        private List<Predicate<Pair<Atom, Atom>>> atomPairPredicates;

        AtomPairSelection(AtomContainer atomContainer) {
            this.atomPairs = Combinatorics.uniquePairsOf(atomContainer.getAtoms())
                    .collect(Collectors.toList());
            this.atomPairPredicates = new ArrayList<>();
        }

        public AtomPairSelection distance(double distanceThreshold) {
            atomPairPredicates.add(pair -> LinearAlgebra3D.distanceFast(pair.getLeft().getCoordinates(),
                    pair.getRight().getCoordinates()) < distanceThreshold * distanceThreshold);
            return this;
        }

        public Stream<Pair<Atom, Atom>> asFilteredAtomPairs() {
            return atomPairs.stream()
                    .filter(Selection.compose(atomPairPredicates));
        }
    }

    public static class GroupPairSelection extends AtomPairSelection {
        private List<Pair<Group, Group>> groupPairs;
        private List<Predicate<Pair<Group, Group>>> groupPairPredicates;

        GroupPairSelection(GroupContainer groupContainer) {
            super(groupContainer);
            this.groupPairs = Combinatorics.uniquePairsOf(groupContainer.getGroups())
                    .collect(Collectors.toList());
            this.groupPairPredicates = new ArrayList<>();
        }

        public GroupPairSelection alphaCarbonDistance(double distanceThreshold) {
            return distance(distanceThreshold, new AlphaCarbonRepresentationScheme());
        }

        public GroupPairSelection betaCarbonDistance(double distanceThreshold) {
            return distance(distanceThreshold, new BetaCarbonRepresentationScheme());
        }

        public GroupPairSelection distance(double distanceThreshold, RepresentationScheme representationScheme) {
            groupPairPredicates.add(pair -> LinearAlgebra3D.distanceFast(representationScheme.determineRepresentingAtom(pair.getLeft()).getCoordinates(),
                    representationScheme.determineRepresentingAtom(pair.getRight()).getCoordinates()) < distanceThreshold * distanceThreshold);
            return this;
        }

        public Stream<Pair<Group, Group>> asFilteredGroupPairs() {
            return groupPairs.stream()
                    .filter(Selection.compose(groupPairPredicates));
        }
    }
}
