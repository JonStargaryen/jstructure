package de.bioforscher.jstructure.feature.topology;

import de.bioforscher.jstructure.feature.FeatureProvider;
import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceAreaCalculator;
import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.model.structure.AminoAcid;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.ResidueContainer;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * <p>An adaptation of the ANVIL algorithm which processes protein structures and places the membrane in the most suitable
 * way.</p>
 *
 * <p>Intellectual and implementation credit to:</p>
 * <p><b>[Postic, 2015] - Membrane positioning for high- and low-resolution protein structures through a binary classification approach - Guillaume Postic, Yassine Ghouzam, Vincent Guiraud, and Jean-Christophe Gelly, Protein Engineering, Design & Selection, 2015, 1–5, doi: 10.1093/protein/gzv063</b></p>
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
public class ANVIL implements FeatureProvider<Protein> {
    /**
     * Definition of amino acids which prefer to interact with water as a solvent and have a low tendency of being
     * embedded in the cell membrane.
     */
    public static final ResidueContainer.AminoAcidFilter SOLVENT_AMINO_ACID_FILTER =
            new ResidueContainer.AminoAcidFilter(Stream.of("ARG", "ASP", "LYS", "GLU", "ASN", "GLN", "PRO", "THR",
                    "TYR").map(AminoAcid::valueOfIgnoreCase)
                    .collect(Collectors.toList()), false);

    private static final int DEFAULT_NUMBER_OF_SPHERE_POINTS = 350;
    private static final double DEFAULT_STEP = 1.0;
    private static final double DEFAULT_MINTHICK = 20.0;
    private static final double DEFAULT_MAXTHICK = 40.0;
    private static final double DEFAULT_AFILTER = 40.0;
    private static final double DEFAULT_DENSITY_OF_MEMBRANE_POINTS = 2.0;

    /**
     * parameters
     */
    private double step;
    private double minthick;
    private double maxthick;
    private double afilter;
    private int numberOfSpherePoints;
    private double density;
    private double maximalExtent;
    private double[] centerOfMass;
    private Protein protein;

    /**
     *
     * @param numberOfSpherePoints
     * @param afilter
     * @param minthick
     * @param maxthick
     * @param step
     * @param density
     */
    public ANVIL(int numberOfSpherePoints, double afilter, double minthick, double maxthick, double step,
                 double density) {
        this.numberOfSpherePoints = numberOfSpherePoints;
        this.afilter = afilter;
        this.minthick = minthick;
        this.maxthick = maxthick;
        this.step = step;
        this.density = density;
    }

    public ANVIL() {
        this(DEFAULT_NUMBER_OF_SPHERE_POINTS, DEFAULT_AFILTER, DEFAULT_MINTHICK, DEFAULT_MAXTHICK, DEFAULT_STEP,
                DEFAULT_DENSITY_OF_MEMBRANE_POINTS);
    }

    @Override
    public void process(Protein featureContainer) {
        this.protein = featureContainer;
        //TODO implement
    }

    /**
     * creates a number of axes close to that of the given membrane
     * @param membrane a well-positioned, but potentially not optimal, membrane
     * @return 350 axes close to that of the membrane
     */
    private List<double[]> findProximateAxes(PotentialMembrane membrane) {
        List<double[]> points = generateSpherePoints(30000);
        Collections.sort(points, new Comparator<double[]>() {
            @Override
            public int compare(double[] d1, double[] d2) {
                return Double.compare(distance(d1), distance(d2));
            }

            private double distance(double[] d) {
                return LinearAlgebra3D.distance(d, membrane.point);
            }
        });
        return points.subList(0, this.numberOfSpherePoints);
    }

    /**
     * counts the total number of exposed hydrophobic respectively hydrophilic residues
     * @return
     */
    private int[] hphobHphil() {
        // delegate the more fine-grained impl when not interested in the placement of a residue relative
        // to the potential membrane plane
        return hphobHphil(false, null, null, null);
    }

    /**
     * counts how well hydrophobic and hydrophilic residues exposed to the solvent are embedded by the membrane
     * @param checkMembranePlane when false residues only need to be exposed in order to count (this is used for the initial 'global' counting of hydrophobic/hydrophilic residues within the structure)
     * @param diam parameters describing the membrane placement
     * @param c1 parameters describing the membrane placement
     * @param c2 parameters describing the membrane placement
     * @return [countOfHydrophobicResidues, countOfHydrophilicResidues]
     */
    private int[] hphobHphil(boolean checkMembranePlane, double[] diam, double[] c1, double[] c2) {
        int[] hphobHphil = { 0, 0 };

        protein.alphaCarbons().filter(residue -> residue.getFeature(double.class,
                    AccessibleSurfaceAreaCalculator.FeatureNames.ACCESSIBLE_SURFACE_AREA.name()) <
                afilter).filter(residue -> true);
        //TODO implement
//        for(Chain chain : this.protein.chains) {
//            for(Residue residue : chain.residues) {
//                // skip residues with too low ASA values - in the original code this is
//                // checked after determining whether the residue is within the membrane
//                // but this check should be way faster and, thus, reduce computational load
//                if(residue.features.get(FeatureType.ACCESSIBLE_SURFACE_AREA.name())[0] < this.afilter) {
//                    continue;
//                }
//
//                // give the option to ignore the membrane placement
//                if(checkMembranePlane && !isInSpace(this.modelConverter.getCA(residue).xyz, diam, c1, c2)) {
//                    continue;
//                }
//
//                if(MEMBRANE_AMINO_ACIDS.contains(residue.aminoAcid)) {
//                    hphobHphil[0]++;
//                } else {
//                    hphobHphil[1]++;
//                }
//            }
//        }

        return hphobHphil;
    }

    /**
     * evaluates the quality of any given membrane slice
     * @param hphil
     * @param hphob
     * @param hphiltotal
     * @param hphobtotal
     * @return the qValue - the higher, the better is the embedding described by this membrane
     */
    private double qValue(double hphil, double hphob, double hphiltotal, double hphobtotal) {
        if(hphobtotal < 1) {
            hphobtotal=0.1;
        }
        if(hphiltotal < 1) {
            hphiltotal += 1;
        }
        double partTotal = hphob + hphil;

        return (hphob * (hphiltotal - hphil) - hphil * (hphobtotal - hphob)) /
                (partTotal * hphobtotal * hphiltotal * (hphobtotal + hphiltotal - partTotal));
    }

    /**
     * return a point S so that S = V * d + P
     * @param vector V
     * @param distance d
     * @param point P
     * @return a point
     */
    private double[] thales(double[] vector, double distance, double[] point) {
        return LinearAlgebra3D.add(point, LinearAlgebra3D.multiply(vector, distance / LinearAlgebra3D.norm(vector)));
    }

    /**
     * generates a defined number of points on a sphere with radiues <code>maximalExtent</code> around <code>centerOfMass</code>
     * @return a collection of sphere points
     */
    private List<double[]> generateSpherePoints(int numberOfSpherePoints) {
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

class PotentialMembrane {
    double[] point;
    double[] c1;
    double[] c2;
    double hphob;
    double hphil;
    double total;
    double qmax;

    public PotentialMembrane(double[] c1, double[] c2, int[] hphobHphil) {
        this.c1 = c1;
        this.c2 = c2;
        this.hphob = hphobHphil[0];
        this.hphil = hphobHphil[1];
        this.total = this.hphob + this.hphil;
    }

    public PotentialMembrane(double[] spherePoint, double[] c1, double[] c2, double hphob, double hphil, double total,
                             double qmax) {
        this.point = spherePoint;
        this.c1 = c1;
        this.c2 = c2;
        this.hphob = hphob;
        this.hphil = hphil;
        this.total = total;
        this.qmax = qmax;
    }
}