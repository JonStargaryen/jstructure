package de.bioforscher.jstructure.feature.topology;

import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceAreaCalculator;
import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.mathematics.LinearAlgebraAtom;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.*;
import java.util.function.IntConsumer;
import java.util.function.Predicate;

/**
 * <p>An adaptation of the ANVIL algorithm which processes protein structures and places the membrane in the most suitable
 * way.</p>
 *
 * <p>Intellectual and implementation credit to:</p>
 * <p><b>[Postic, 2015] - Membrane positioning for high- and low-resolution protein structures through a binary classification approach - Guillaume Postic, Yassine Ghouzam, Vincent Guiraud, and Jean-Christophe Gelly, Protein Engineering, Design & SelectionOld, 2015, 1–5, doi: 10.1093/protein/gzv063</b></p>
 *
 * <pre>© Univ. Paris Diderot & Inserm, 2015
 guillaume.postic@univ-paris-diderot.fr
 This software is a computer program whose purpose is to assign membrane
 boundaries to a protein three-dimensional structure, by using the spatial
 coordinates of the alpha carbons. It can also process coarse-grained
 models of protein structures. The algorithm follows an approach that
 treats the problem of membrane assignment as a binary classification.
 The software is implemented in Python and requires the PyPy interpreter.
 It also requires Naccess 2.1.1 (S.J. Hubbard, 1996) to calculate the
 atomic accessible surface.
 This software is governed by the CeCILL license under French law and
 abiding by the rules of distribution of free software.  You can  use,
 modify and/ or redistribute the software under the terms of the CeCILL
 license as circulated by CEA, CNRS and INRIA at the following URL
 "http://www.cecill.info".
 As a counterpart to the access to the source code and  rights to copy,
 modify and redistribute granted by the license, users are provided only
 with a limited warranty  and the software's author,  the holder of the
 economic rights,  and the successive licensors  have only  limited
 liability.
 In this respect, the user's attention is drawn to the risks associated
 with loading,  using,  modifying and/or developing or reproducing the
 software by the user in light of its specific status of free software,
 that may mean  that it is complicated to manipulate,  and  that  also
 therefore means  that it is reserved for developers  and  experienced
 professionals having in-depth computer knowledge. Users are therefore
 encouraged to load and test the software's suitability as regards their
 requirements in conditions enabling the security of their systems and/or
 data to be ensured and,  more generally, to use and operate it in the
 same conditions as regards security.
 The fact that you are presently reading this means that you have had
 knowledge of the CeCILL license and that you accept its terms.</pre>
 */
@FeatureProvider(provides = { ANVIL.MEMBRANE, ANVIL.TOPOLOGY }, requires = { AccessibleSurfaceAreaCalculator.ACCESSIBLE_SURFACE_AREA })
public class ANVIL extends AbstractFeatureProvider {
    private static final Logger logger = LoggerFactory.getLogger(ANVIL.class);
    public static final String MEMBRANE = "MEMBRANE";
    public static final String TOPOLOGY = "TOPOLOGY";

    /*
     * default values
     */
    private static final int DEFAULT_NUMBER_OF_SPHERE_POINTS = 350;
    private static final double DEFAULT_STEP = 1.0;
    private static final double ALTERNATIVE_STEP = 0.3;
    private static final double DEFAULT_MINTHICK = 20.0;
    private static final double DEFAULT_MAXTHICK = 40.0;
    private static final double DEFAULT_AFILTER = 40.0;
    private static final double DEFAULT_DENSITY_OF_MEMBRANE_POINTS = 2.0;

    /*
     * parameters
     */
    private double minthick;
    private double maxthick;
    private int numberOfSpherePoints;
    private double density;
    private Predicate<Group> asaFilter;

    /**
     * The fine-grained constructor.
     * @param numberOfSpherePoints how many points to generate initially
     * @param afilter the maximal accessible surface area value of a getResidue to still be considered
     * @param minthick the minimum thickness of the membrane
     * @param maxthick the maximum thickness of the membrane
     * @param density how dense pseudo atoms representing the membrane layers are placed
     */
    public ANVIL(int numberOfSpherePoints, double afilter, double minthick, double maxthick, double density) {
        this.numberOfSpherePoints = numberOfSpherePoints;
        this.minthick = minthick;
        this.maxthick = maxthick;
        this.density = density;
        this.asaFilter = residue -> residue.getFeatureAsDouble(AccessibleSurfaceAreaCalculator.ACCESSIBLE_SURFACE_AREA) > afilter;
    }

    /**
     * The constructor using default values.
     */
    public ANVIL() {
        this(DEFAULT_NUMBER_OF_SPHERE_POINTS,
                DEFAULT_AFILTER,
                DEFAULT_MINTHICK,
                DEFAULT_MAXTHICK,
                DEFAULT_DENSITY_OF_MEMBRANE_POINTS);
    }

    @Override
    protected void processInternally(Protein protein) {
        //TODO generally bugged

        AtomContainer alphaCarbons = Selection.on(protein)
                .alphaCarbonAtoms()
                .asAtomContainer();
        // compute center of mass based on alpha carbons which equals centroid() since every atoms weights the same
        double[] centerOfMass = LinearAlgebraAtom.centroid(alphaCarbons);
        double maximalExtent = 1.2 * LinearAlgebraAtom.maximalExtent(alphaCarbons, centerOfMass);

        int[] initialHphobHphil = hphobHphil(protein);
        Membrane initialMembrane = processSpherePoints(protein,
                generateSpherePoints(numberOfSpherePoints, centerOfMass, maximalExtent),
                centerOfMass,
                DEFAULT_STEP,
                initialHphobHphil);
        logger.debug("initial membrane quality: {}", initialMembrane.qmax);

        Membrane alternativeMembrane = processSpherePoints(protein,
                findProximateAxes(initialMembrane, centerOfMass, maximalExtent),
                centerOfMass,
                DEFAULT_STEP,
                initialHphobHphil);
        logger.debug("alternative membrane quality: {}", alternativeMembrane.qmax);

        Membrane membrane = initialMembrane.qmax > alternativeMembrane.qmax ? initialMembrane : alternativeMembrane;

        logger.debug("membrane thickness is {} A", LinearAlgebra3D.distance(membrane.planePoint1, membrane.planePoint2));
        //TODO adjustment of thickness
        logger.debug("embedding quality is {}", membrane.qmax);

        assignTopology(protein, membrane);
        placeMembraneMolecules(protein, membrane, maximalExtent);
    }

    /**
     * Places evenly distributed pseudo-atoms for sake of membrane visualization.
     */
    private void placeMembraneMolecules(Protein protein, Membrane membrane, double maximalExtent) {
        // square maximal extent to speed up computations later on
        double radius = maximalExtent * maximalExtent;
        double[] normalVector = membrane.normalVector;
        for(double[] layer : Arrays.asList(membrane.planePoint1, membrane.planePoint2)) {
            double d = - LinearAlgebra3D.dotProduct(normalVector, layer);
            for(double i = -1000; i < 1000; i += density) {
                for(double j = -1000; j < 1000; j += density) {
                    double[] atom = new double[] { i, j, 0 };
                    try {
                        atom[2] = -(d + i * normalVector[0] + j * normalVector[1]) / normalVector[2];
                    } catch (ArithmeticException e) {
                        logger.warn(e.getLocalizedMessage());
                    }

                    // distance cutoff is also squared
                    if(LinearAlgebra3D.distanceFast(atom, layer) <= radius &&
                            minimalSquaredDistanceToProteinCAAtom(protein, atom) > 12.0) {
                        membrane.membraneAtoms.add(atom);
                    }
                }
            }
        }
    }

    /**
     * Computes the minimal distance of a point to any alpha carbon of a structure.
     * @param membranePseudoAtom the point to check
     * @return the closest distance to any alpha carbon of the protein
     */
    private double minimalSquaredDistanceToProteinCAAtom(Protein protein, final double[] membranePseudoAtom) {
        return Selection.on(protein)
                    .aminoAcids()
                    .alphaCarbonAtoms()
                    .asFilteredAtoms()
                    .map(Atom::getCoordinates)
                    .mapToDouble(coordinates -> LinearAlgebra3D.distanceFast(coordinates, membranePseudoAtom))
                    .min()
                    .orElse(Double.MAX_VALUE);
    }

    /**
     * Assign the topology to each getResidue which is embedded in the membrane.
     */
    private void assignTopology(Protein protein, Membrane membrane) {
        protein.setFeature(MEMBRANE, membrane);
        Selection.on(protein)
                .aminoAcids()
                .asFilteredGroups()
                .filter(residue -> isInMembranePlane(Selection.on(residue)
                                .alphaCarbonAtoms()
                                .asAtom()
                                .getCoordinates(),
                        membrane.normalVector,
                        membrane.planePoint1,
                        membrane.planePoint2))
                .forEach(residue ->  residue.setFeature(TOPOLOGY, membrane)
        );
    }

    /**
     * Checks a set of points to check and returns the most suitable membrane orientation.
     * @param spherePoints a collection of sphere points to evaluate
     * @return the most promising membrane for the given set of points
     */
    private Membrane processSpherePoints(Protein protein, List<double[]> spherePoints, double[] centerOfMass, double step, int[] initialHphobHphil) {
        // best performing membrane
        Membrane membrane = null;
        // best performing membrane's score
        double qmax = 0;

        // construct slices of thickness 1.0 along the axis connecting the centerOfMass and the spherePoint
        for (double[] spherePoint : spherePoints) {
            double[] diam = LinearAlgebra3D.multiply(LinearAlgebra3D.subtract(centerOfMass, spherePoint), 2.0);
            double diamNorm = LinearAlgebra3D.norm(diam);

            List<Membrane> qvartemp = new ArrayList<>();

            for (double i = 0; i < diamNorm - step; i += step) {
                double dPointC1 = i;
                double dPointC2 = i + step;

                double[] c1 = thales(diam, dPointC1, spherePoint);
                double[] c2 = thales(diam, dPointC2, spherePoint);

                // evaluate how well this membrane slice embedded the peculiar residues
                int[] hphobHphil = hphobHphil(protein, true, diam, c1, c2);

                qvartemp.add(new Membrane(c1, c2, hphobHphil[0], hphobHphil[1]));
            }

            int jmax = (int) ((minthick / step) - 1);

            for (double width = 0; width < maxthick; width = (jmax + 1) * step) {
                int imax = qvartemp.size() - 1 - jmax;

                for (int i = 0; i < imax; i++) {
                    double[] c1 = qvartemp.get(i).planePoint1;
                    double[] c2 = qvartemp.get(i + jmax).planePoint2;

                    double hphob = 0;
                    double hphil = 0;
                    double total = 0;

                    for (int j = 0; j < jmax; j++) {
                        Membrane ij = qvartemp.get(i + j);
                        if (j == 0 || j == jmax - 1) {
                            hphob += 0.5 * ij.hphob;
                            hphil += 0.5 * ij.hphil;
                        } else {
                            hphob += ij.hphob;
                            hphil += ij.hphil;
                        }
                        total += ij.total;
                    }

                    if (hphob > 0) {
                        double qvaltest = qValue(hphil, hphob, initialHphobHphil[1], initialHphobHphil[0]);
                        if (qvaltest > qmax) {
                            qmax = qvaltest;
                            membrane = new Membrane(spherePoint, c1, c2, hphob, hphil, total, qmax);
                        }
                    }
                }
                jmax++;
            }
        }

        return membrane;
    }

    /**
     * Creates a number of axes close to that of the given membrane.
     * @param membrane a well-positioned, but potentially not optimal, membrane
     * @return 30000 axes close to that of the membrane
     */
    private List<double[]> findProximateAxes(Membrane membrane, double[] centerOfMass, double maximalExtent) {
        List<double[]> points = generateSpherePoints(30000, centerOfMass, maximalExtent);
        points.sort(new Comparator<double[]>() {
            @Override
            public int compare(double[] d1, double[] d2) {
                return Double.compare(distanceFast(d1), distanceFast(d2));
            }

            private double distanceFast(double[] d) {
                return LinearAlgebra3D.distanceFast(d, membrane.point);
            }
        });
        return points.subList(0, this.numberOfSpherePoints);
    }

    /**
     * Counts the total number of exposed hydrophobic respectively hydrophilic residues.
     * @return an <code>int[]</code>: the first index is the count of hydrophobic, the second that of hydrophilic
     */
    private int[] hphobHphil(Protein protein) {
        // delegate the more fine-grained impl when not interested in the placement of a getResidue relative
        // to the potential membrane plane
        return hphobHphil(protein, false, null, null, null);
    }

    /**
     * counts how well hydrophobic and hydrophilic residues exposed to the solvent are embedded by the membrane
     * @param checkMembranePlane when false residues only need to be exposed in order to count (this is used for the initial 'global' counting of hydrophobic/hydrophilic residues within the structure)
     * @param diam parameters describing the membrane placement
     * @param c1 parameters describing the membrane placement
     * @param c2 parameters describing the membrane placement
     * @return [countOfHydrophobicResidues, countOfHydrophilicResidues]
     */
    private int[] hphobHphil(Protein protein, final boolean checkMembranePlane, final double[] diam, final double[] c1, final double[] c2) {
        return Selection.on(protein)
                .aminoAcids()
                .asFilteredGroups()
                // select out exposed residues
                .filter(asaFilter)
                // select for residues within the membrane plane
                .filter(residue -> !checkMembranePlane ||
                      isInMembranePlane(Selection.on(residue)
                              .alphaCarbonAtoms()
                              .asAtom()
                              .getCoordinates(), diam, c1, c2))
                // map to index representation - 0: membrane-tendency, 1: polar getResidue
                .mapToInt(residue -> residue.getGroupInformation().getAminoAcidFamily().getANVILGrouping().ordinal())
                // unknown amino acids could be mapped to indices exceeding the array
                .filter(index -> index < 2)
                .collect(HphobHphilConsumer::new, HphobHphilConsumer::accept, HphobHphilConsumer::combine)
                .getHphobHpil();
    }

    static class HphobHphilConsumer implements IntConsumer {
        private int[] values = new int[2];

        int[] getHphobHpil() {
            return values;
        }

        public void accept(int index) {
            values[index]++;
        }

        void combine(HphobHphilConsumer other) {
            values[0] += other.values[0];
            values[1] += other.values[1];
        }
    }

    /**
     * Checks whether a certain point is within the membrane planes.
     * @param pointToTest the peculiar points
     * @param normalVector the normal vector of this membrane
     * @param planePoint1 a point representing one side of the membrane
     * @param planePoint2 the point representing the other side of the membrane
     * @return true iff the given point is within the membrane
     */
    private boolean isInMembranePlane(double[] pointToTest, double[] normalVector, double[] planePoint1, double[] planePoint2) {
        normalVector = LinearAlgebra3D.normalize(normalVector);
        final double d1 = - LinearAlgebra3D.dotProduct(normalVector, planePoint1);
        final double d2 = - LinearAlgebra3D.dotProduct(normalVector, planePoint2);
        final double d = - LinearAlgebra3D.dotProduct(normalVector, pointToTest);
        return d > Math.min(d1, d2) && d < Math.max(d1, d2);
    }

    /**
     * Evaluates the quality of any given membrane slice.
     * @param hphil the current number of embedded hydrophilic residues
     * @param hphob the current number of embedded hydrophobic residues
     * @param hphiltotal the maximal number of hydrophilic residues
     * @param hphobtotal the maximal number of hydrophobic residues
     * @return the qValue - the higher, the better is the embedding described by this membrane
     */
    private double qValue(double hphil, double hphob, double hphiltotal, double hphobtotal) {
        if(hphobtotal < 1) {
            hphobtotal = 0.1;
        }
        if(hphiltotal < 1) {
            hphiltotal += 1;
        }
        double partTotal = hphob + hphil;

        return (hphob * (hphiltotal - hphil) - hphil * (hphobtotal - hphob)) /
                (partTotal * hphobtotal * hphiltotal * (hphobtotal + hphiltotal - partTotal));
    }

    /**
     * Return a point S so that S = V * d + P.
     * @param vector V
     * @param distance d
     * @param point P
     * @return a point
     */
    private double[] thales(double[] vector, double distance, double[] point) {
        return LinearAlgebra3D.add(point, LinearAlgebra3D.multiply(vector, distance / LinearAlgebra3D.norm(vector)));
    }

    /**
     * Generates a defined number of points on a sphere with radiues <code>maximalExtent</code> around
     * <code>centerOfMass</code>
     * @return a collection of sphere points
     */
    private List<double[]> generateSpherePoints(int numberOfSpherePoints, double[] centerOfMass, double maximalExtent) {
        List<double[]> points = new ArrayList<>();

        double oldPhi = 0;
        for(int k = 1; k < numberOfSpherePoints + 1; k++) {
            double h = -1 + 2 * (k - 1) / ((double) (numberOfSpherePoints - 1));
            double theta = Math.acos(h);
            double phi = (k == 1 || k == numberOfSpherePoints) ? 0 : (oldPhi + 3.6 / Math.sqrt(numberOfSpherePoints*(1 -
                    h * h))) % (2 * Math.PI);

            double[] point = new double[] {
                    maximalExtent * Math.sin(phi) * Math.sin(theta) + centerOfMass[0],
                    maximalExtent * Math.cos(theta) + centerOfMass[1],
                    maximalExtent * Math.cos(phi) * Math.sin(theta) + centerOfMass[2]
            };
            points.add(point);
            oldPhi = phi;
        }

        return points;
    }
}
