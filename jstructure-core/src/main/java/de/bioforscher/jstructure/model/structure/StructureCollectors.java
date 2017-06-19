package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.mathematics.LinearAlgebra;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;
import java.util.stream.Collector;

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

    public static Collector<Atom, ?, BasicAtomContainer> toAtomContainer() {
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

        BasicAtomContainer getContainer() {
            return new BasicAtomContainer(atoms,
                    atoms.get(0).getParentGroup().getParentChain().getParentProtein());
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

    public static Collector<Group, ?, BasicGroupContainer> toGroupContainer() {
        return Collector.of(GroupContainerConsumer::new,
                GroupContainerConsumer::accept,
                GroupContainerConsumer::combine,
                GroupContainerConsumer::getContainer,
                Collector.Characteristics.CONCURRENT);
    }

    public static class GroupContainerConsumer implements Consumer<Group> {
        List<Group> groups;

        GroupContainerConsumer() {
            this.groups = new ArrayList<>();
        }

        BasicGroupContainer getContainer() {
            return new BasicGroupContainer(groups,
                    groups.get(0).getParentChain().getParentProtein());
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

    public static Collector<Chain, ?, BasicChainContainer> toChainContainer() {
        return Collector.of(ChainContainerConsumer::new,
                ChainContainerConsumer::accept,
                ChainContainerConsumer::combine,
                ChainContainerConsumer::getContainer,
                Collector.Characteristics.CONCURRENT);
    }

    public static class ChainContainerConsumer implements Consumer<Chain> {
        List<Chain> chains;

        ChainContainerConsumer() {
            this.chains = new ArrayList<>();
        }

        BasicChainContainer getContainer() {
            return new BasicChainContainer(chains,
                    chains.get(0).getParentProtein());
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
