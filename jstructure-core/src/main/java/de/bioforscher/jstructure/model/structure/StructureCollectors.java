package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.mathematics.LinearAlgebra;
import de.bioforscher.jstructure.model.identifier.ProteinIdentifier;
import de.bioforscher.jstructure.model.structure.container.StructureContainer;

import java.util.*;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.stream.Collector;
import java.util.stream.Collectors;

/**
 * A custom collection of structure collectors.
 * Created by S on 09.11.2016.
 */
public class StructureCollectors {
    public static Collector<Atom, ?, double[]> toCentroid() {
        return Collector.of(CentroidAverage::new,
                CentroidAverage::accept,
                CentroidAverage::combine,
                CentroidAverage::average,
                Collector.Characteristics.CONCURRENT);
    }

    static class CentroidAverage implements Consumer<Atom> {
        private double[] total;
        private int count;

        CentroidAverage() {
            this.total = new double[3];
            this.count = 0;
        }

        double[] average() {
            return count > 0 ? LinearAlgebra.on(total).divide(count).getValue() : new double[3];
        }

        @Override
        public void accept(Atom atom) {
            total = LinearAlgebra.on(total).add(atom.getCoordinates()).getValue();
            count++;
        }

        CentroidAverage combine(CentroidAverage other) {
            total = LinearAlgebra.on(total).add(other.total).getValue();
            count += other.count;
            return this;
        }
    }

    public static Collector<Atom, ?, double[]> toCenterOfMass() {
        return Collector.of(CenterOfMassAverage::new,
                CenterOfMassAverage::accept,
                CenterOfMassAverage::combine,
                CenterOfMassAverage::average,
                Collector.Characteristics.CONCURRENT);
    }

    static class CenterOfMassAverage implements Consumer<Atom> {
        private double[] coordinate;
        private double mass;

        CenterOfMassAverage() {
            this.coordinate = new double[3];
            this.mass = 0;
        }

        double[] average() {
            return mass > 0 ? LinearAlgebra.on(coordinate).divide(mass).getValue() : new double[3];
        }

        @Override
        public void accept(Atom atom) {
            coordinate = LinearAlgebra.on(atom.getCoordinates()).multiply(atom.getElement().getAtomicMass()).add(coordinate).getValue();
            mass += atom.getElement().getAtomicMass();
        }

        CenterOfMassAverage combine(CenterOfMassAverage other) {
            coordinate = LinearAlgebra.on(coordinate).add(other.coordinate).getValue();
            mass += other.mass;
            return this;
        }
    }

    public static <S extends StructureContainer> Collector<S, ?, Structure> toIsolatedStructure() {
        return Collector.of(IsolatedStructureConsumer::new,
                IsolatedStructureConsumer::accept,
                IsolatedStructureConsumer::combine,
                IsolatedStructureConsumer::getIsolatedStructure);
    }

    static class IsolatedStructureConsumer<S extends StructureContainer> implements Consumer<S> {
        List<S> elements;

        IsolatedStructureConsumer() {
            this.elements = new ArrayList<>();
        }

        @SuppressWarnings("unchecked")
        Structure getIsolatedStructure() {
            // return empty container when nothing was selected
            if(elements.isEmpty()) {
                return new Structure(ProteinIdentifier.UNKNOWN_PROTEIN_IDENTIFIER);
            }

            // determine the lists content based on the first element, this can crash horribly when the list is for some
            // reason populated with elements of varying type
            S referenceElement = elements.get(0);
            // if element is whole structure, perform mere deep copy
            if(referenceElement instanceof Structure) {
                return (Structure) referenceElement.createDeepCopy();
            }

            // non-trivial cases: basic procedure is to deep copy all children and link them to the a shallow copy of
            // their parent

            // the case for chains
            if(referenceElement instanceof Chain) {
                Chain referenceChain = (Chain) referenceElement;
                Structure clonedStructure = referenceChain.getParentStructure().createShallowCopy();

                ((List<Chain>) elements)
                        .stream()
                        .map(Chain::createDeepCopy)
                        .forEach(clonedStructure::addChain);

                return clonedStructure;
            }

            if(referenceElement instanceof Group) {
                Group referenceGroup = (Group) referenceElement;
                Structure clonedStructure = referenceGroup.getParentChain().getParentStructure().createShallowCopy();

                // keep track of association
                Map<Chain, List<Group>> parentChainAssociation = ((List<Group>) elements).stream()
                        .collect(Collectors.groupingBy(Group::getParentChain,
                                LinkedHashMap::new,
                                Collectors.mapping(Function.identity(), Collectors.toList())));
                for(Map.Entry<Chain, List<Group>> entry : parentChainAssociation.entrySet()) {
                    Chain clonedChain = entry.getKey().createShallowCopy();
                    List<Group> clonedGroups = entry.getValue().stream()
                            .map(Group::createDeepCopy)
                            .collect(Collectors.toList());
                    // link current chain to the cloned structure
                    clonedStructure.addChain(clonedChain);
                    // every deep-cloned group points to the shallow copy of the current chain
                    clonedGroups.forEach(clonedChain::addGroup);
                }

                return clonedStructure;
            }
            if(referenceElement instanceof Atom) {
                Atom referenceAtom = (Atom) referenceElement;
                Structure clonedStructure = referenceAtom.getParentGroup().getParentChain().getParentStructure().createShallowCopy();

                // keep track of association
                Map<Group, List<Atom>> parentGroupAssociation = ((List<Atom>) elements).stream()
                        .collect(Collectors.groupingBy(Atom::getParentGroup,
                                LinkedHashMap::new,
                                Collectors.mapping(Function.identity(), Collectors.toList())));
                // key: original chain, value: shallow clone of chain
                Map<Chain, Chain> clonedParentChains = new HashMap<>();
                for (Map.Entry<Group, List<Atom>> entry : parentGroupAssociation.entrySet()) {
                    Group originalGroup = entry.getKey();
                    Group clonedGroup = originalGroup.createShallowCopy();
                    List<Atom> clonedAtoms = entry.getValue().stream()
                            .map(Atom::createDeepCopy)
                            .collect(Collectors.toList());

                    // if map does not yet contain clone of the parent chain, create it
                    Chain originalParentChain = originalGroup.getParentChain();
                    clonedParentChains.putIfAbsent(originalParentChain, originalParentChain.createShallowCopy());
                    Chain clonedChain = clonedParentChains.get(originalParentChain);

                    // link current group to the cloned chain
                    clonedChain.addGroup(clonedGroup);
                    // every deep-cloned atom points to the shallow copy of the current group
                    clonedAtoms.forEach(clonedGroup::addAtom);
                }

                // link all chain copies to the cloned structure
                clonedParentChains.values()
                        .forEach(clonedStructure::addChain);

                return clonedStructure;
            }

            throw new IllegalArgumentException("could not create isolated container as parent references were null");
        }

        @Override
        public void accept(S element) {
            elements.add(element);
        }

        IsolatedStructureConsumer<S> combine(IsolatedStructureConsumer<S> other) {
            elements.addAll(other.elements);
            return this;
        }
    }

//    public static <S extends StructureContainer> Collector<S, ?, Structure> toSubStructure() {
//        return Collector.of(SubStructureConsumer::new,
//                SubStructureConsumer::accept,
//                SubStructureConsumer::combine,
//                SubStructureConsumer::getSubStructure,
//                Collector.Characteristics.CONCURRENT);
//    }
//
//    static class SubStructureConsumer<S extends StructureContainer> implements Consumer<S> {
//        List<S> elements;
//
//        SubStructureConsumer() {
//            this.elements = new ArrayList<>();
//        }
//
//        @SuppressWarnings("unchecked")
//        Structure getSubStructure() {
//            // return empty container when nothing was selected
//            if(elements.isEmpty()) {
//                return new Structure(ProteinIdentifier.UNKNOWN_PROTEIN_IDENTIFIER);
//            }
//
//            // determine the lists content based on the first element, this can crash horribly when the list is for some
//            // reason populated with elements of varying type
//            S referenceElement = elements.get(0);
//            // if element is whole structure, return it
//            if(referenceElement instanceof Structure) {
//                return (Structure) referenceElement;
//            }
//
//            // non-trivial cases: basic procedure is to create a representation of the structure containing only
//            // selected elements
//
//            // the case for chains
//            if(referenceElement instanceof Chain) {
//                Chain referenceChain = (Chain) referenceElement;
//                Structure clonedStructure = referenceChain.getParentStructure().createShallowCopy();
//
//                ((List<Chain>) elements)
//                        .stream()
//                        .map(Chain::createDeepCopy)
//                        .forEach(clonedStructure::addChain);
//
//                return clonedStructure;
//            }
//
//            if(referenceElement instanceof Group) {
//                Group referenceGroup = (Group) referenceElement;
//                Structure clonedStructure = referenceGroup.getParentChain().getParentStructure().createShallowCopy();
//
//                // keep track of association
//                Map<Chain, List<Group>> parentChainAssociation = ((List<Group>) elements).stream()
//                        .collect(Collectors.groupingBy(Group::getParentChain));
//                for(Map.Entry<Chain, List<Group>> entry : parentChainAssociation.entrySet()) {
//                    Chain clonedChain = entry.getKey().createShallowCopy();
//                    List<Group> clonedGroups = entry.getValue().stream()
//                            .map(Group::createDeepCopy)
//                            .collect(Collectors.toList());
//                    // link current chain to the cloned structure
//                    clonedStructure.addChain(clonedChain);
//                    // every deep-cloned group points to the shallow copy of the current chain
//                    clonedGroups.forEach(clonedChain::addGroup);
//                }
//
//                return clonedStructure;
//            }
//            if(referenceElement instanceof Atom) {
//                Atom referenceAtom = (Atom) referenceElement;
//                Structure clonedStructure = referenceAtom.getParentGroup().getParentChain().getParentStructure().createShallowCopy();
//
//                // keep track of association
//                Map<Group, List<Atom>> parentGroupAssociation = ((List<Atom>) elements).stream()
//                        .collect(Collectors.groupingBy(Atom::getParentGroup));
//                // key: original chain, value: shallow clone of chain
//                Map<Chain, Chain> clonedParentChains = new HashMap<>();
//                for (Map.Entry<Group, List<Atom>> entry : parentGroupAssociation.entrySet()) {
//                    Group originalGroup = entry.getKey();
//                    Group clonedGroup = originalGroup.createShallowCopy();
//                    List<Atom> clonedAtoms = entry.getValue().stream()
//                            .map(Atom::createDeepCopy)
//                            .collect(Collectors.toList());
//
//                    // if map does not yet contain clone of the parent chain, create it
//                    Chain originalParentChain = clonedGroup.getParentChain();
//                    clonedParentChains.putIfAbsent(originalParentChain, originalParentChain.createShallowCopy());
//                    Chain clonedChain = clonedParentChains.get(originalParentChain);
//
//                    // link current group to the cloned chain
//                    clonedChain.addGroup(clonedGroup);
//                    // every deep-cloned atom points to the shallow copy of the current group
//                    clonedAtoms.forEach(clonedGroup::addAtom);
//                }
//
//                // link all chain copies to the cloned structure
//                clonedParentChains.values()
//                        .forEach(clonedStructure::addChain);
//
//                return clonedStructure;
//            }
//
//            throw new IllegalArgumentException("could not create isolated container as parent references were null");
//        }
//
//        @Override
//        public void accept(S element) {
//            elements.add(element);
//        }
//
//        SubStructureConsumer<S> combine(SubStructureConsumer<S> other) {
//            elements.addAll(other.elements);
//            return this;
//        }
//    }
}
