package de.bioforscher.jstructure.model;

import de.bioforscher.jstructure.alignment.AlignmentAlgorithm;
import de.bioforscher.jstructure.alignment.AlignmentResult;
import de.bioforscher.jstructure.alignment.SVDSuperimposer;
import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.container.ChainContainer;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.filter.AtomNameFilter;
import org.apache.commons.math3.exception.DimensionMismatchException;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;
import java.util.stream.Collector;
import java.util.stream.Collectors;

/**
 * A custom collection of structure collectors.
 * Created by S on 09.11.2016.
 */
public class StructureCollectors {
    public static Collector<Protein, ?, List<Protein>> toAlignedEnsemble() {
        return Collector.of(ProteinSuperimposerByReference::new,
                ProteinSuperimposerByReference::accept,
                ProteinSuperimposerByReference::combine,
                ProteinSuperimposerByReference::getAlignedProteins,
                Collector.Characteristics.CONCURRENT);
    }

    static class ProteinSuperimposerByReference implements Consumer<Protein> {
        /**
         * the protein all entries are aligned against
         */
        Protein reference;
        /**
         * the container of aligned proteins
         */
        List<Protein> alignedProteins;
        /**
         * the approach to align fragments
         */
        AlignmentAlgorithm alignmentStrategy;

        ProteinSuperimposerByReference() {
            alignedProteins = new ArrayList<>();
            alignmentStrategy = new SVDSuperimposer(AtomNameFilter.BACKBONE_ATOM_FILTER);
        }

        @Override
        public void accept(Protein protein) {
            if(reference == null) {
                reference = protein;
            }

            alignedProteins.add(protein);
            try {
                AlignmentResult alignmentResult = alignmentStrategy.align(reference, protein);
                alignmentResult.transform(protein);
            } catch (DimensionMismatchException e) {
                e.printStackTrace();
            }
        }

        ProteinSuperimposerByReference combine(ProteinSuperimposerByReference other) {
            other.alignedProteins.forEach(this);
            return this;
        }

        List<Protein> getAlignedProteins() {
            return alignedProteins;
        }
    }

    public static Collector<Atom, ?, double[]> toCentroid() {
        return Collector.of(CentroidAverager::new,
                CentroidAverager::accept,
                CentroidAverager::combine,
                CentroidAverager::average,
                Collector.Characteristics.CONCURRENT);
    }

    static class CentroidAverager implements Consumer<Atom> {
        private double[] total;
        private int count;

        CentroidAverager() {
            this.total = new double[3];
            this.count = 0;
        }

        double[] average() {
            return count > 0 ? LinearAlgebra3D.divide(total, count) : new double[3];
        }

        @Override
        public void accept(Atom atom) {
            total = LinearAlgebra3D.add(total, atom.getCoordinates());
            count++;
        }

        CentroidAverager combine(CentroidAverager other) {
            total = LinearAlgebra3D.add(total, other.total);
            count += other.count;
            return this;
        }
    }

    public static Collector<Atom, ?, double[]> toCenterOfMass() {
        return Collector.of(CenterOfMassAverager::new,
                CenterOfMassAverager::accept,
                CenterOfMassAverager::combine,
                CenterOfMassAverager::average,
                Collector.Characteristics.CONCURRENT);
    }

    static class CenterOfMassAverager implements Consumer<Atom> {
        private double[] coordinate;
        private double mass;

        CenterOfMassAverager() {
            this.coordinate = new double[3];
            this.mass = 0;
        }

        double[] average() {
            return mass > 0 ? LinearAlgebra3D.divide(coordinate, mass) : new double[3];
        }

        @Override
        public void accept(Atom atom) {
            coordinate = LinearAlgebra3D.add(coordinate, LinearAlgebra3D.multiply(atom.getCoordinates(),
                    atom.getElement().getAtomicMass()));
            mass += atom.getElement().getAtomicMass();
        }

        CenterOfMassAverager combine(CenterOfMassAverager other) {
            coordinate = LinearAlgebra3D.add(coordinate, other.coordinate);
            mass += other.mass;
            return this;
        }
    }

    public static Collector<Atom, ?, AtomContainer> toAtomContainer() {
        return Collector.of(AtomContainerConsumer::new,
                AtomContainerConsumer::accept,
                AtomContainerConsumer::combine,
                AtomContainerConsumer::getContainer,
                Collector.Characteristics.CONCURRENT);
    }

    static class AtomContainerConsumer implements Consumer<Atom> {
        List<Atom> atoms;

        AtomContainerConsumer() {
            this.atoms = new ArrayList<>();
        }

        AtomContainer getContainer() {
            return new BasicAtomContainer(atoms);
        }

        @Override
        public void accept(Atom atom) {
            atoms.add(atom);
        }

        AtomContainerConsumer combine(AtomContainerConsumer other) {
            atoms.addAll(other.atoms);
            return this;
        }
    }

    static class BasicAtomContainer implements AtomContainer {
        private List<Atom> atoms;
        private Map<Enum, Object> featureMap;

        BasicAtomContainer(List<Atom> atomList) {
            this.atoms = atomList;
            this.featureMap = new HashMap<>();
        }

        @Override
        public List<Atom> getAtoms() {
            return atoms;
        }

        @Override
        public Map<Enum, Object> getFeatureMap() {
            return featureMap;
        }
    }

    public static Collector<Group, ?, GroupContainer> toGroupContainer() {
        return Collector.of(GroupContainerConsumer::new,
                GroupContainerConsumer::accept,
                GroupContainerConsumer::combine,
                GroupContainerConsumer::getContainer,
                Collector.Characteristics.CONCURRENT);
    }

    static class GroupContainerConsumer implements Consumer<Group> {
        List<Group> groups;

        GroupContainerConsumer() {
            this.groups = new ArrayList<>();
        }

        GroupContainer getContainer() {
            return new BasicGroupContainer(groups);
        }

        @Override
        public void accept(Group group) {
            groups.add(group);
        }

        GroupContainerConsumer combine(GroupContainerConsumer other) {
            groups.addAll(other.groups);
            return this;
        }
    }

    static class BasicGroupContainer implements AtomContainer, GroupContainer {
        private List<Group> groups;
        private Map<Enum, Object> featureMap;

        BasicGroupContainer(List<Group> groupList) {
            this.groups = groupList;
            this.featureMap = new HashMap<>();
        }

        @Override
        public List<Atom> getAtoms() {
            return groups().flatMap(Group::atoms).collect(Collectors.toList());
        }

        @Override
        public List<Group> getGroups() {
            return groups;
        }

        @Override
        public Map<Enum, Object> getFeatureMap() {
            return featureMap;
        }
    }

    public static Collector<Chain, ?, ChainContainer> toChainContainer() {
        return Collector.of(ChainContainerConsumer::new,
                ChainContainerConsumer::accept,
                ChainContainerConsumer::combine,
                ChainContainerConsumer::getContainer,
                Collector.Characteristics.CONCURRENT);
    }

    static class ChainContainerConsumer implements Consumer<Chain> {
        List<Chain> chains;

        ChainContainerConsumer() {
            this.chains = new ArrayList<>();
        }

        ChainContainer getContainer() {
            return new BasicChainContainer(chains);
        }

        @Override
        public void accept(Chain chain) {
            chains.add(chain);
        }

        ChainContainerConsumer combine(ChainContainerConsumer other) {
            chains.addAll(other.chains);
            return this;
        }
    }

    static class BasicChainContainer implements ChainContainer {
        private List<Chain> chains;
        private Map<Enum, Object> featureMap;

        BasicChainContainer(List<Chain> atomList) {
            this.chains = atomList;
            this.featureMap = new HashMap<>();
        }

        @Override
        public List<Chain> getChains() {
            return chains;
        }

        @Override
        public List<Group> getGroups() {
            return chains().flatMap(Chain::groups).collect(Collectors.toList());
        }

        @Override
        public List<Atom> getAtoms() {
            return groups().flatMap(Group::atoms).collect(Collectors.toList());
        }

        @Override
        public Map<Enum, Object> getFeatureMap() {
            return featureMap;
        }
    }
}
