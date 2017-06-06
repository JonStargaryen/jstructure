package de.bioforscher.jstructure.feature.topology;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.structure.Group;

import java.util.Collection;
import java.util.List;

/**
 * Represents the membrane layer of proteins as a set of 2 layers of points.
 * Created by bittrich on 5/17/17.
 */
public class MembraneContainer extends FeatureContainerEntry {
    private List<double[]> membraneAtoms;
    private List<TransMembraneSubunit> transMembraneHelices;
    private String hydrophobicThickness;
    private String tiltAngle;
    private String deltaGTransfer;
    private String topology;

    MembraneContainer(AbstractFeatureProvider featureProvider) {
        super(featureProvider);
    }

    public void addTransMembraneHelix(TransMembraneSubunit transMembraneSubunit) {
        transMembraneHelices.add(transMembraneSubunit);
    }

    public List<double[]> getMembraneAtoms() {
        return membraneAtoms;
    }

    public void setMembraneAtoms(List<double[]> membraneAtoms) {
        this.membraneAtoms = membraneAtoms;
    }

    public List<TransMembraneSubunit> getTransMembraneHelices() {
        return transMembraneHelices;
    }

    public void setTransMembraneHelices(List<TransMembraneSubunit> transMembraneHelices) {
        this.transMembraneHelices = transMembraneHelices;
    }

    public String getHydrophobicThickness() {
        return hydrophobicThickness;
    }

    public void setHydrophobicThickness(String hydrophobicThickness) {
        this.hydrophobicThickness = hydrophobicThickness;
    }

    public String getTiltAngle() {
        return tiltAngle;
    }

    public void setTiltAngle(String tiltAngle) {
        this.tiltAngle = tiltAngle;
    }

    public String getDeltaGTransfer() {
        return deltaGTransfer;
    }

    public void setDeltaGTransfer(String deltaGTransfer) {
        this.deltaGTransfer = deltaGTransfer;
    }

    public String getTopology() {
        return topology;
    }

    public void setTopology(String topology) {
        this.topology = topology;
    }

    public boolean isTransmembraneGroup(Group group) {
        int residueNumber = group.getResidueNumber().getResidueNumber();
        return transMembraneHelices.stream()
                .filter(transMembraneSubunit -> transMembraneSubunit.getChainId().equals(group.getParentChain().getChainId().getChainId()))
                .map(TransMembraneSubunit::getSegments)
                .flatMap(Collection::stream)
                .anyMatch(range -> range.getLeft() <= residueNumber && range.getRight() >= residueNumber);
    }

    //TODO a class for the topology of individual groups
}
