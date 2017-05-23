package de.bioforscher.jstructure.feature.topology;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;

import java.util.List;

/**
 * Represents the membrane layer of proteins as a set of 2 layers of points.
 * Created by bittrich on 5/17/17.
 */
public class Membrane extends FeatureContainerEntry {
    private List<double[]> membraneAtoms;
    private List<TransMembraneHelix> transMembraneHelices;
    private String hydrophobicThickness;
    private String tiltAngle;
    private String deltaGTransfer;
    private String topology;

    Membrane(AbstractFeatureProvider featureProvider) {
        super(featureProvider);
    }

    public void addTransMembraneHelix(TransMembraneHelix transMembraneHelix) {
        transMembraneHelices.add(transMembraneHelix);
    }

    public List<double[]> getMembraneAtoms() {
        return membraneAtoms;
    }

    public void setMembraneAtoms(List<double[]> membraneAtoms) {
        this.membraneAtoms = membraneAtoms;
    }

    public List<TransMembraneHelix> getTransMembraneHelices() {
        return transMembraneHelices;
    }

    public void setTransMembraneHelices(List<TransMembraneHelix> transMembraneHelices) {
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

    //TODO checks if some residue is trans-membrane or not
}
