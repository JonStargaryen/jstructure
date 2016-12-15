package de.bioforscher.jstructure.feature.topology;

import de.bioforscher.jstructure.feature.motif.SequenceMotif;
import de.bioforscher.jstructure.mathematics.CoordinateManipulations;
import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.model.structure.AminoAcid;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.*;
import java.util.function.Function;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * The container for membranes.
 * Created by S on 03.11.2016.
 */
public class Membrane {
    private static final String AMINO_ACID_POTENTIAL = "topology/potential/amino_acid_potential.dat";
    /**
     * the collection of points representing the membrane layers
     */
    final List<double[]> membraneAtoms;
    /**
     * the map of functions assigning a potential to each amino acid observation
     */
    private static Map<AminoAcid, Function<Double, Double>> fittingFunctions;
    /**
     * the center point of the membrane
     */
    double[] point;
    /**
     * the border point on one side of the membrane
     */
    double[] planePoint1;
    /**
     * the border point on other side of the membrane
     */
    double[] planePoint2;
    double hphob;
    double hphil;
    double total;
    /**
     * the quality of the embedding described by this membrane
     */
    double qmax;
    /**
     * the normal vector of this membrane
     */
    double[] normalVector;

    public Membrane(double[] c1, double[] c2, double hphob, double hphil) {
        this.planePoint1 = c1;
        this.planePoint2 = c2;
        this.hphob = hphob;
        this.hphil = hphil;
        this.total = hphob + hphil;
        this.normalVector = LinearAlgebra3D.subtract(c1, c2);
        this.point = LinearAlgebra3D.subtract(planePoint1, LinearAlgebra3D.multiply(normalVector, 0.5));
        this.membraneAtoms = new ArrayList<>();
        if(fittingFunctions == null) {
            initializeLibrary();
        }
    }

    public Membrane(double[] point, double[] c1, double[] c2, double hphob, double hphil, double total, double qmax) {
        this(c1, c2, hphob, hphil);
        this.point = point;
        this.total = total;
        this.qmax = qmax;
    }

    private InputStream getResourceAsStream(String filepath) {
        return Objects.requireNonNull(Thread.currentThread().getContextClassLoader().getResourceAsStream(filepath),
                "failed to findAny resource as InputStream");
    }

    private synchronized void initializeLibrary() {
        InputStream idxIs = getResourceAsStream(AMINO_ACID_POTENTIAL);
        fittingFunctions = new BufferedReader(new InputStreamReader(idxIs)).lines()
                // skip intial header line
                .filter(line -> !line.startsWith("Residue"))
                .map(FittingFunction::new)
                .collect(Collectors.toMap(FittingFunction::getAminoAcid, Function.identity()));
    }

    public double getQmax() {
        return qmax;
    }

    /**
     * Assigns a statically derived energy term to each amino acid based on its distance to the center of the membrane.
     *
     * see reference:
     * Proteins. 2005 May 1;59(2):252-65.
     * Properties of integral membrane protein structures: derivation of an implicit membrane potential.
     * Ulmschneider MB, Sansom MS, Di Nola A.
     */
    static class FittingFunction implements Function<Double, Double> {
        AminoAcid aminoAcid;
        double[] parameters;

        static final FittingFunction DEFAULT_FITTING_FUNCTION = new FittingFunction("UNK 0 0 0 0 0 0 0");

        FittingFunction(String line) {
            this.aminoAcid = AminoAcid.valueOfIgnoreCase(line.substring(0, 3));
            this.parameters = Pattern.compile("\\s+").splitAsStream(line)
                    // skip amino acid annotation
                    .skip(1)
                    // skip extra values at the end
                    .limit(7)
                    .mapToDouble(Double::valueOf)
                    .toArray();
        }

        @Override
        public Double apply(Double distanceToCenterOfMembrane) {
            double q1 = distanceToCenterOfMembrane - parameters[3];
            double q2 = distanceToCenterOfMembrane - parameters[6];
            return parameters[0] +
                    parameters[1] * Math.exp(-parameters[2] * q1 * q1) +
                    parameters[4] * Math.exp(-parameters[5] * q2 * q2);
        }

        AminoAcid getAminoAcid() {
            return aminoAcid;
        }
    }

    /**
     * Computes the potential for a sequence motif. Each individual potential is computed and averaged.
     * @param sequenceMotif a collection of residues
     * @return average energy observed for this
     */
    public double computePotential(SequenceMotif sequenceMotif) {
       return sequenceMotif
               .getGroupContainer()
               .groups()
               .mapToDouble(this::computePotential)
               .average()
               .getAsDouble();
    }

    /**
     * Computes the individual potential.
     * @param residue a collection of residues
     * @return observed potential
     */
    public double computePotential(Group residue) {
        return fittingFunctions.getOrDefault(AminoAcid.valueOfIgnoreCase(residue.getPdbName()),
                FittingFunction.DEFAULT_FITTING_FUNCTION).apply(distanceToMembraneCenter(residue));
    }

    /**
     * Computes the distance of a given point to the plane trough the center of the described membrane.
     * @param coordinates some point
     * @return the distance to the center of the membrane plane
     */
    public double distanceToMembraneCenter(double[] coordinates) {
        double[] w = LinearAlgebra3D.normalize(normalVector);
        return LinearAlgebra3D.dotProduct(coordinates, w) + LinearAlgebra3D.dotProduct(point, w);
    }

    public double distanceToMembraneCenter(Atom atom) {
        return distanceToMembraneCenter(atom.getCoordinates());
    }

    /**
     * Computes the distance of a given point to the plane through the center of the described membrane. The collection
     * of atoms is represented by its centroid.
     * @param atomContainer some collection of atoms
     * @return the distance to the center of te membrane plane
     */
    public double distanceToMembraneCenter(AtomContainer atomContainer) {
        return distanceToMembraneCenter(CoordinateManipulations.centroid(atomContainer));
    }
}