package de.bioforscher.jstructure.model.structure.selection;

import de.bioforscher.jstructure.mathematics.IntegerInterval;
import de.bioforscher.jstructure.mathematics.LinearAlgebra;
import de.bioforscher.jstructure.mathematics.Pair;
import de.bioforscher.jstructure.model.identifier.ChainIdentifier;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
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

        void inferContainerName(AtomContainer container) {
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

            return prefilteredAtomStream
                    .filter(Selection.compose(atomPredicates));
        }

        /**
         * @see StructureCollectors#toIsolatedStructure() for documentation of this call's behaviour
         */
        public Structure asIsolatedStructure() {
            Structure structure = asFilteredAtoms()
                    .collect(StructureCollectors.toIsolatedStructure());
            inferContainerName(structure);
            return structure;
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

        public AtomSelection pdbSerial(IntegerInterval... pdbSerialRanges) {
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

            return prefilteredGroupStream
                    .filter(Selection.compose(groupPredicates));
        }

        public Optional<Group> asOptionalGroup() {
            return asFilteredGroups()
                    .findFirst();
        }

        public Optional<AminoAcid> asOptionalAminoAcid() {
            Optional<Group> optional = asFilteredGroups()
                    .findFirst();
            try {
                return optional
                        .map(AminoAcid.class::cast);
            } catch (ClassCastException e) {
                // wrap potential ClassCastException as SelectionException
                throw new SelectionException(optional + " cannot be cast to " + AminoAcid.class.getSimpleName(), e);
            }
        }

        public Group asGroup() {
            return asOptionalGroup()
                    .orElseThrow(() -> new SelectionException("did not find group matching " + Stream.of(groupPredicates,
                            chainPredicates)
                            .flatMap(Collection::stream)
                            .map(Pair::getRight)
                            .collect(Collectors.joining(", ", "[",  "]")) + " in " + groupContainer.getIdentifier()));
        }

        public AminoAcid asAminoAcid() {
            return asOptionalAminoAcid()
                    .orElseThrow(() -> new SelectionException("did not find amino acid matching " + Stream.of(groupPredicates,
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
                    .contains(group.getResidueIdentifier().getResidueNumber()), "residue numbers: " + Arrays.toString(residueNumbers));
            return this;
        }

        public GroupSelection residueNumber(IntegerInterval... residueNumberRanges) {
            registerGroupPredicate(group -> Stream.of(residueNumberRanges)
                    .anyMatch(range -> group.getResidueIdentifier().getResidueNumber() >= range.getLeft() &&
                            group.getResidueIdentifier().getResidueNumber() <= range.getRight()),
                    "residue ranges: " + Arrays.toString(residueNumberRanges));
            return this;
        }

        public GroupSelection residueIdentifier(ResidueIdentifier residueIdentifier) {
            registerGroupPredicate(group -> group.getResidueIdentifier().equals(residueIdentifier),
                    "residueIdentifier: " + residueIdentifier);
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
            return chainContainer.chains()
                    .filter(Selection.compose(chainPredicates));
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

        public ChainSelection chainName(String... chainIds) {
            registerChainPredicate(chain -> Stream.of(chainIds)
                    .anyMatch(chainName -> chainName.equals(chain.getChainIdentifier().getChainId())),
                    "chain names: " + Arrays.toString(chainIds));
            return this;
        }

        public ChainSelection chainId(String chainId) {
            registerChainPredicate(chain -> chain.getChainIdentifier().getChainId().equals(chainId),
                    "chainId: " + chainId);
            return this;
        }

        public ChainSelection chainIdentifier(ChainIdentifier chainIdentifier, ChainIdentifier furtherChainIdentifiers) {
            Set<String> chainIds = Stream.of(chainIdentifier, furtherChainIdentifiers)
                    .map(ChainIdentifier::getChainId)
                    .collect(Collectors.toSet());
            registerChainPredicate(chain -> chainIds.contains(chain.getChainIdentifier().getChainId()),
                    "chainIds: " + chainIds);
            return this;
        }

        public ChainSelection chainIdentifier(ChainIdentifier chainIdentifier) {
            return chainId(chainIdentifier.getChainId());
        }
    }
}
