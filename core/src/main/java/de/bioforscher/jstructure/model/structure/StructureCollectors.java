package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.container.ChainContainer;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;
import java.util.stream.Collector;

/**
 * A custom collection of structure collectors.
 * Created by S on 09.11.2016.
 */
public class StructureCollectors {
    private static final Logger logger = LoggerFactory.getLogger(StructureCollectors.class);

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
            return new Group(atoms);
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
            return new Chain(groups);
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
            return new Protein(chains);
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
}
