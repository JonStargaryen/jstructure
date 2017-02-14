package de.bioforscher.jstructure.model.structure.selection;

import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.model.structure.*;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.container.ChainContainer;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.container.StructureContainer;
import de.bioforscher.jstructure.model.structure.family.AminoAcidFamily;

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

    public static GroupSelection on(GroupContainer groupContainer) {
        return new GroupSelection(Objects.requireNonNull(groupContainer));
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
        String specifiedName;
        String parentContainerName;
        boolean negationMode;

        AtomSelection(AtomContainer atomContainer) {
            this.atomPredicates = new ArrayList<>();
            this.groupPredicates = new ArrayList<>();
            this.chainPredicates = new ArrayList<>();
            this.atomContainer = atomContainer;
            this.parentContainerName = atomContainer.getIdentifier();
        }

        /**
         * Request this selection to clone all returned elements so they can be manipulated independently of their
         * original reference.
         * @return the builder
         */
        public AtomSelection cloneElements() {
            cloneRequested = true;
            return this;
        }

        /**
         * Request this selection to negate all coming predicates.
         * @return the builder
         */
        public AtomSelection negationModeEnter() {
            this.negationMode = true;
            return this;
        }

        /**
         * Request this selection to negate all coming predicates.
         * @return the builder
         */
        public AtomSelection negationModeLeave() {
            this.negationMode = false;
            return this;
        }

        /**
         * Allows to specifiedName the returned selection. I.e. {@link StructureContainer#getIdentifier()} will return this
         * functions argument and will not try to infer the specifiedName itself from the given query.
         * @param name the specifiedName this selection shall have
         * @return the builder
         */
        public AtomSelection nameContainer(String name) {
            this.specifiedName = name;
            return this;
        }

        void inferContainerName(StructureContainer container) {
            if(specifiedName != null) {
                container.setIdentifier(specifiedName);
                return;
            }

            container.setIdentifier(parentContainerName);
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
            AtomContainer container = asFilteredAtoms()
                    .collect(StructureCollectors.toAtomContainer());

            inferContainerName(container);

            return container;
        }

        public Optional<Atom> asOptionalAtom() {
            return asFilteredAtoms()
                    .findAny();
        }

        public Atom asAtom() {
            return asOptionalAtom()
                    .orElseThrow(NoSuchElementException::new);
        }

        private void registerAtomPredicate(Predicate<Atom> atomPredicate) {
            atomPredicates.add(negationMode ? atomPredicate.negate() : atomPredicate);
        }

        public AtomSelection element(Collection<Element> elements) {
            return element(elements.toArray(new Element[elements.size()]));
        }

        public AtomSelection element(Element... elements) {
            registerAtomPredicate(atom -> Stream.of(elements)
                    .anyMatch(element -> atom.getElement().equals(element)));
            return this;
        }

        public AtomSelection atomName(Collection<String> atomNames){
            return atomName(atomNames.toArray(new String[atomNames.size()]));
        }

        public AtomSelection atomName(String... atomNames) {
            registerAtomPredicate(atom -> Stream.of(atomNames)
                    .anyMatch(atomName -> atomName.equals(atom.getName())));
            return this;
        }

        public AtomSelection backboneAtoms() {
            registerAtomPredicate(atom -> AminoAcidFamily.ATOM_NAMES.BACKBONE_ATOM_NAMES.contains(atom.getName()));
            return this;
        }

        public AtomSelection sideChainAtoms() {
            registerAtomPredicate(atom -> AminoAcidFamily.ATOM_NAMES.SIDECHAIN_ATOM_NAMES.contains(atom.getName()));
            return this;
        }

        public AtomSelection allAtoms() {
            registerAtomPredicate(atom -> Stream.concat(AminoAcidFamily.ATOM_NAMES.BACKBONE_ATOM_NAMES.stream(),
                    AminoAcidFamily.ATOM_NAMES.SIDECHAIN_ATOM_NAMES.stream()).collect(Collectors.toList()).contains(atom.getName()));
            return this;
        }

        public AtomSelection alphaCarbonAtoms() {
            registerAtomPredicate(atom -> AminoAcidFamily.ATOM_NAMES.CA_ATOM_NAME.equals(atom.getName()));
            return this;
        }

        public AtomSelection backboneCarbonAtoms() {
            registerAtomPredicate(atom -> AminoAcidFamily.ATOM_NAMES.C_ATOM_NAME.equals(atom.getName()));
            return this;
        }

        public AtomSelection backboneNitrogenAtoms() {
            registerAtomPredicate(atom -> AminoAcidFamily.ATOM_NAMES.N_ATOM_NAME.equals(atom.getName()));
            return this;
        }

        public AtomSelection backboneOxygenAtoms() {
            registerAtomPredicate(atom -> AminoAcidFamily.ATOM_NAMES.O_ATOM_NAME.equals(atom.getName()));
            return this;
        }

        public AtomSelection betaCarbonAtoms() {
            registerAtomPredicate(atom -> AminoAcidFamily.ATOM_NAMES.CB_ATOM_NAME.equals(atom.getName()));
            return this;
        }

        public AtomSelection hydrogenAtoms() {
            registerAtomPredicate(atom -> AminoAcidFamily.ATOM_NAMES.H_ATOM_NAMES.contains(atom.getName()));
            return this;
        }

        public AtomSelection backboneHydrogen() {
//            backboneAtoms(); - 18/1/17 fix: will never be true as backbone atoms no longer container hydrogen names
            hydrogenAtoms();
            return this;
        }

        public AtomSelection nonHydrogenAtoms() {
            registerAtomPredicate(atom -> !AminoAcidFamily.ATOM_NAMES.H_ATOM_NAMES.contains(atom.getName()));
            return this;
        }

        public AtomSelection pdbSerial(int... pdbSerials) {
            registerAtomPredicate(atom -> IntStream.of(pdbSerials).boxed().collect(Collectors.toList())
                    .contains(atom.getPdbSerial()));
            return this;
        }

        public AtomSelection pdbSerial(IntegerRange... pdbSerialRanges) {
            registerAtomPredicate(atom -> Stream.of(pdbSerialRanges)
                    .anyMatch(range -> range.getLeft() >= atom.getPdbSerial() && range.getRight() <= atom.getPdbSerial()));
            return this;
        }

        public AtomSelection distance(Atom referenceAtom, double distanceCutoff) {
            return distance(referenceAtom.getCoordinates(), distanceCutoff);
        }

        public AtomSelection distance(double[] coordinates, double distanceCutoff) {
            registerAtomPredicate(atom -> LinearAlgebra3D.distanceFast(coordinates,
                    atom.getCoordinates()) < distanceCutoff * distanceCutoff);
            return this;
        }
    }

    public static class GroupSelection extends AtomSelection {
        GroupSelection(GroupContainer groupContainer) {
            super(groupContainer);
            this.groupContainer = groupContainer;
            this.parentContainerName = groupContainer.getIdentifier();
        }


        public AtomSelection atomSelection() {
            return this;
        }

        @Override
        public GroupSelection cloneElements() {
            super.cloneElements();
            return this;
        }

        @Override
        public GroupSelection nameContainer(String name) {
            super.nameContainer(name);
            return this;
        }

        @Override
        public GroupSelection negationModeEnter() {
            super.negationModeEnter();
            return this;
        }

        @Override
        public GroupSelection negationModeLeave() {
            super.negationModeLeave();
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
            GroupContainer container = asFilteredGroups()
                    .collect(StructureCollectors.toGroupContainer());

            inferContainerName(container);

            return container;
        }

        public Optional<Group> asOptionalGroup() {
            return asFilteredGroups()
                    .findAny();
        }

        public Group asGroup() {
            return asOptionalGroup()
                    .orElseThrow(NoSuchElementException::new);
        }

        private void registerGroupPredicate(Predicate<Group> groupPredicate) {
            groupPredicates.add(negationMode ? groupPredicate.negate() : groupPredicate);
        }

        public GroupSelection aminoAcids() {
            registerGroupPredicate(Group::isAminoAcid);
            return this;
        }

        public GroupSelection aminoAcids(AminoAcidFamily... aminoAcids) {
            // pre-filter for group type amino acid
            aminoAcids();
            registerGroupPredicate(group -> Stream.of(aminoAcids)
                    .map(AminoAcidFamily::getThreeLetterCode)
                    .anyMatch(threeLetterCode -> threeLetterCode.equals(group.getThreeLetterCode())));
            return this;
        }

        public GroupSelection nucleotides() {
            registerGroupPredicate(Group::isNucleotide);
            return this;
        }

        public GroupSelection hetatms() {
            registerGroupPredicate(Group::isLigand);
            return this;
        }

        public GroupSelection water() {
            registerGroupPredicate(group -> group.getThreeLetterCode().equals("HOH"));
            return this;
        }

        public GroupSelection groupName(String... groupNames) {
            registerGroupPredicate(group -> Stream.of(groupNames)
                    .anyMatch(groupName -> groupName.equals(group.getThreeLetterCode())));
            return this;
        }

        public GroupSelection residueNumber(int... residueNumbers) {
            registerGroupPredicate(group -> IntStream.of(residueNumbers).boxed().collect(Collectors.toList())
                    .contains(group.getResidueNumber()));
            return this;
        }

        public GroupSelection residueNumber(IntegerRange... residueNumberRanges) {
            registerGroupPredicate(group -> Stream.of(residueNumberRanges)
                    .anyMatch(range -> group.getResidueNumber() >= range.getLeft() && group.getResidueNumber() <= range.getRight()));
            return this;
        }

        @Override
        public GroupSelection distance(Atom referenceAtom, double distanceCutoff) {
            return distance(referenceAtom.getCoordinates(), distanceCutoff);
        }

        @Override
        public GroupSelection distance(double[] coordinates, double distanceCutoff) {
            registerGroupPredicate(group -> group.atoms()
                    .anyMatch(atom -> LinearAlgebra3D.distanceFast(coordinates, atom.getCoordinates()) < distanceCutoff * distanceCutoff));
            return this;
        }
    }

    public static class ChainSelection extends GroupSelection {
        ChainSelection(ChainContainer chainContainer) {
            super(chainContainer);
            this.chainContainer = chainContainer;
            this.parentContainerName = chainContainer.getIdentifier();
        }

        public GroupSelection groupSelection() {
            return this;
        }

        @Override
        public ChainSelection cloneElements() {
            super.cloneElements();
            return this;
        }

        @Override
        public ChainSelection nameContainer(String name) {
            super.nameContainer(name);
            return this;
        }

        @Override
        public ChainSelection negationModeEnter() {
            super.negationModeEnter();
            return this;
        }

        @Override
        public ChainSelection negationModeLeave() {
            super.negationModeLeave();
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
            ChainContainer container = asFilteredChains()
                    .collect(StructureCollectors.toChainContainer());

            inferContainerName(container);

            return container;
        }

        public Optional<Chain> asOptionalChain() {
            return asFilteredChains()
                    .findAny();
        }

        public Chain asChain() {
            return asOptionalChain()
                    .orElseThrow(NoSuchElementException::new);
        }

        private void registerChainPredicate(Predicate<Chain> chainPredicate) {
            chainPredicates.add(negationMode ? chainPredicate.negate() : chainPredicate);
        }

        public ChainSelection chainName(String... chainNames) {
            registerChainPredicate(chain -> Stream.of(chainNames)
                    .anyMatch(chainName -> chainName.equals(chain.getChainId())));
            return this;
        }
    }
}
