package de.bioforscher.jstructure.model.structure.selection;

import de.bioforscher.jstructure.mathematics.LinearAlgebra;
import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.structure.*;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.container.ChainContainer;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.container.StructureContainer;

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
    private static <T> Predicate<T> compose(List<Pair<Predicate<T>, String>> predicates) {
        return predicates.stream()
                // the actual predicate is on the left
                .map(Pair::getLeft)
                .reduce(Predicate::and)
                .orElse(x -> true);
    }

    public static class AtomSelection {
        AtomContainer atomContainer;
        GroupContainer groupContainer;
        ChainContainer chainContainer;
        List<Pair<Predicate<Atom>, String>> atomPredicates;
        List<Pair<Predicate<Group>, String>> groupPredicates;
        List<Pair<Predicate<Chain>, String>> chainPredicates;
        static final String CUSTOM_PREDICATE = "custom predicate";
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

        public AtomSelection customAtomPredicate(Predicate<Atom> atomPredicate) {
            registerAtomPredicate(atomPredicate, CUSTOM_PREDICATE);
            return this;
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
                    .findFirst();
        }

        public Atom asAtom() {
            return asOptionalAtom()
                    .orElseThrow(() -> new SelectionException("did not find atom matching " + Stream.of(atomPredicates,
                            groupPredicates,
                            chainPredicates)
                            .flatMap(Collection::stream)
                            .map(Pair::getRight)
                            .collect(Collectors.joining(", ", "[",  "]")) + " in " + atomContainer.getIdentifier()));
        }

        private void registerAtomPredicate(Predicate<Atom> atomPredicate, String description) {
            atomPredicates.add(new Pair<>(negationMode ? atomPredicate.negate() : atomPredicate,
                    (negationMode ? "NOT: " : "") + description));
        }

        public AtomSelection element(Collection<Element> elements) {
            return element(elements.toArray(new Element[elements.size()]));
        }

        public AtomSelection element(Element... elements) {
            registerAtomPredicate(atom -> Stream.of(elements)
                    .anyMatch(element -> atom.getElement().equals(element)), "element of: " + Arrays.toString(elements));
            return this;
        }

        public AtomSelection atomName(Collection<String> atomNames){
            return atomName(atomNames.toArray(new String[atomNames.size()]));
        }

        public AtomSelection atomName(String... atomNames) {
            registerAtomPredicate(atom -> Stream.of(atomNames)
                    .anyMatch(atomName -> atomName.equals(atom.getName())), "atom name of: " + Arrays.toString(atomNames));
            return this;
        }

        public AtomSelection alphaCarbonAtoms() {
            registerAtomPredicate(atom -> AminoAcid.ALPHA_CARBON_NAME.equals(atom.getName()), "alpha carbons");
            return this;
        }

        public AtomSelection backboneCarbonAtoms() {
            registerAtomPredicate(atom -> AminoAcid.BACKBONE_CARBON_NAME.equals(atom.getName()), "backbone carbons");
            return this;
        }

        public AtomSelection backboneNitrogenAtoms() {
            registerAtomPredicate(atom -> AminoAcid.BACKBONE_NITROGEN_NAME.equals(atom.getName()), "backbone nitrogens");
            return this;
        }

        public AtomSelection backboneOxygenAtoms() {
            registerAtomPredicate(atom -> AminoAcid.BACKBONE_OXYGEN_NAME.equals(atom.getName()), "backbone oxygens");
            return this;
        }

        public AtomSelection betaCarbonAtoms() {
            registerAtomPredicate(atom -> AminoAcid.BETA_CARBON_NAME.equals(atom.getName()), "beta carbons");
            return this;
        }

        public AtomSelection hydrogenAtoms() {
            registerAtomPredicate(atom -> Group.HYDROGEN_NAMES.contains(atom.getName()), "hydrogens");
            return this;
        }

        public AtomSelection backboneHydrogen() {
//            backboneAtoms(); - 18/1/17 fix: will never be true as backbone atoms no longer container hydrogen names
            hydrogenAtoms();
            return this;
        }

        public AtomSelection nonHydrogenAtoms() {
            registerAtomPredicate(atom -> !Group.HYDROGEN_NAMES.contains(atom.getName()), "non-hydrogens");
            return this;
        }

        public AtomSelection pdbSerial(int... pdbSerials) {
            registerAtomPredicate(atom -> IntStream.of(pdbSerials).boxed().collect(Collectors.toList())
                    .contains(atom.getPdbSerial()), "pdb-serials: " + Arrays.toString(pdbSerials));
            return this;
        }

        public AtomSelection pdbSerial(IntegerRange... pdbSerialRanges) {
            registerAtomPredicate(atom -> Stream.of(pdbSerialRanges)
                    .anyMatch(range -> range.getLeft() >= atom.getPdbSerial() && range.getRight() <= atom.getPdbSerial()),
                    "pdb-serial ranges: " + Arrays.toString(pdbSerialRanges));
            return this;
        }

        public AtomSelection atomDistance(Atom referenceAtom, double distanceCutoff) {
            return atomDistance(referenceAtom.getCoordinates(), distanceCutoff);
        }

        public AtomSelection atomDistance(double[] coordinates, double distanceCutoff) {
            double squaredDistanceCutoff = distanceCutoff * distanceCutoff;
            registerAtomPredicate(atom -> LinearAlgebra.on(coordinates).distanceFast(atom.getCoordinates()) <
                    squaredDistanceCutoff, distanceCutoff + " A around " + Arrays.toString(coordinates));
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

        public GroupSelection customGroupPredicate(Predicate<Group> groupPredicate) {
            registerGroupPredicate(groupPredicate, CUSTOM_PREDICATE);
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
                prefilteredGroupStream = chainContainer.chains()
                        .filter(Selection.compose(chainPredicates))
                        .flatMap(Chain::groups);
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
                    .findFirst();
        }

        public Group asGroup() {
            return asOptionalGroup()
                    .orElseThrow(() -> new SelectionException("did not find group matching " + Stream.of(groupPredicates,
                            chainPredicates)
                            .flatMap(Collection::stream)
                            .map(Pair::getRight)
                            .collect(Collectors.joining(", ", "[",  "]")) + " in " + groupContainer.getIdentifier()));
        }

        private void registerGroupPredicate(Predicate<Group> groupPredicate, String description) {
            groupPredicates.add(new Pair<>(negationMode ? groupPredicate.negate() : groupPredicate, description));
        }

        /**
         * All amino acids which are part of the protein chain (i.e. ligands are ignored).
         * @return a stream of all entities which are part of the protein chain
         */
        public GroupSelection aminoAcids() {
            registerGroupPredicate(Group::isAminoAcid, "amino acids");
            return this;
        }

        public GroupSelection nucleotides() {
            registerGroupPredicate(Group::isNucleotide, "nucleotides");
            return this;
        }

        public GroupSelection hetatms() {
            registerGroupPredicate(Group::isLigand, "ligands");
            return this;
        }

        public GroupSelection water() {
            registerGroupPredicate(Group::isWater, "waters");
            return this;
        }

        public GroupSelection ligands() {
            registerGroupPredicate(Group::isLigand, "ligands");
            return this;
        }

        public GroupSelection groupName(String... groupNames) {
            registerGroupPredicate(group -> Stream.of(groupNames)
                    .anyMatch(groupName -> groupName.equals(group.getThreeLetterCode())), "group names: " + Arrays.toString(groupNames));
            return this;
        }

        public GroupSelection residueNumber(int... residueNumbers) {
            registerGroupPredicate(group -> IntStream.of(residueNumbers).boxed().collect(Collectors.toList())
                    .contains(group.getResidueNumber().getResidueNumber()), "residue numbers: " + Arrays.toString(residueNumbers));
            return this;
        }

        //TODO support for insertion codes

        public GroupSelection residueNumber(IntegerRange... residueNumberRanges) {
            registerGroupPredicate(group -> Stream.of(residueNumberRanges)
                    .anyMatch(range -> group.getResidueNumber().getResidueNumber() >= range.getLeft() &&
                            group.getResidueNumber().getResidueNumber() <= range.getRight()),
                    "residue ranges: " + Arrays.toString(residueNumberRanges));
            return this;
        }

        public GroupSelection groupDistance(Atom referenceAtom, double distanceCutoff) {
            return groupDistance(referenceAtom.getCoordinates(), distanceCutoff);
        }

        public GroupSelection groupDistance(double[] coordinates, double distanceCutoff) {
            double squaredDistanceCutoff = distanceCutoff * distanceCutoff;
            registerGroupPredicate(group -> group.atoms()
                    .anyMatch(atom -> LinearAlgebra.on(coordinates).distanceFast(atom.getCoordinates()) <
                            squaredDistanceCutoff), distanceCutoff + " A around " + Arrays.toString(coordinates));
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

        public ChainSelection customChainPredicate(Predicate<Chain> chainPredicate) {
            registerChainPredicate(chainPredicate, CUSTOM_PREDICATE);
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
                    .findFirst();
        }

        public Chain asChain() {
            return asOptionalChain()
                    .orElseThrow(() -> new SelectionException("did not find chain matching " + Stream.of(chainPredicates)
                            .flatMap(Collection::stream)
                            .map(Pair::getRight)
                            .collect(Collectors.joining(", ", "[",  "]"))  + " in " + chainContainer.getIdentifier()));
        }

        private void registerChainPredicate(Predicate<Chain> chainPredicate, String description) {
            chainPredicates.add(new Pair<>(negationMode ? chainPredicate.negate() : chainPredicate, description));
        }

        public ChainSelection chainName(String... chainNames) {
            registerChainPredicate(chain -> Stream.of(chainNames)
                    .anyMatch(chainName -> chainName.equals(chain.getChainId().getChainId())), "chain names: " + Arrays.toString(chainNames));
            return this;
        }
    }
}
