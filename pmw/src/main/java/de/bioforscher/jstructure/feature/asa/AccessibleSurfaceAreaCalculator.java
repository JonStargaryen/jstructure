package de.bioforscher.jstructure.feature.asa;

import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Element;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * <p>Computes the accessible surface area of each getResidue in a {@link Protein}.</p>
 *
 * <p>Wide parts are BioJava code. Original BioJava doc:</p>
 * <pre>Class to calculate Accessible Surface Areas based on
 * the rolling ball algorithm by Shrake and Rupley.
 *
 * The code is adapted from a python implementation at http://boscoh.com/protein/asapy
 * (now source is available at https://github.com/boscoh/asa).
 * Thanks to Bosco K. Ho for a great piece of code and for his fantastic blog.
 *
 * See
 * Shrake, A., and J. A. Rupley. "Environment and Exposure to Solvent of Protein Atoms.
 * Lysozyme and Insulin." JMB (1973) 79:351-371.
 * Lee, B., and Richards, F.M. "The interpretation of Protein Structures: Estimation of
 * Static Accessibility" JMB (1971) 55:379-400</pre>
 * @author duarte_j
 */
public class AccessibleSurfaceAreaCalculator implements FeatureProvider {
    private static final Logger logger = LoggerFactory.getLogger(AccessibleSurfaceAreaCalculator.class);
    // Bosco uses as default 960, Shrake and Rupley seem to use in their paper 92 (not sure if this is actually the same
    // parameter)
    public static final int DEFAULT_N_SPHERE_POINTS = 960;
    public static final double DEFAULT_PROBE_SIZE = 1.4;
    // Chothia's amino acid atoms vdw radii
    public static final double TRIGONAL_CARBON_VDW = 1.76;
    public static final double TETRAHEDRAL_CARBON_VDW = 1.87;
    public static final double TRIGONAL_NITROGEN_VDW = 1.65;
    public static final double TETRAHEDRAL_NITROGEN_VDW = 1.50;
    // values diverge from the vdw in enum Element
    public static final double SULFUR_VDW = 1.85;
    public static final double OXIGEN_VDW = 1.40;

    public enum FeatureNames {
        ATOM_RADIUS,
        ACCESSIBLE_SURFACE_AREA
    }

    private int numberOfSpherePoints;
    private double probeSize;
    private double cons;
    private AtomContainer nonHydrogenAtoms;
    private List<double[]> spherePoints;

    /**
     * The fine-grained constructor.
     * @param numberOfSpherePoints how many points to use
     * @param probeSize sphere size to probe for ASA
     */
    public AccessibleSurfaceAreaCalculator(int numberOfSpherePoints, double probeSize) {
        this.numberOfSpherePoints = numberOfSpherePoints;
        this.probeSize = probeSize;
    }

    /**
     * The default constructor using default parameters.
     */
    public AccessibleSurfaceAreaCalculator() {
        this(DEFAULT_N_SPHERE_POINTS, DEFAULT_PROBE_SIZE);
    }

    private void assignRadius(Atom atom) {
        atom.setFeature(FeatureNames.ATOM_RADIUS, determineRadius(atom));
    }

    private void assignSingleAsa(Group group) {
        group.setFeature(FeatureNames.ACCESSIBLE_SURFACE_AREA,
                Selection.on(group)
                        .nonHydrogenAtoms()
                        .asFilteredAtoms()
                        .mapToDouble(this::calcSingleAsa)
                        .sum());
    }

    @Override
    public void process(Protein protein) {
        // determine radius for all non-hydrogen atoms and assign it to the atom's internal feature map
        GroupContainer residueContainer = Selection.on(protein)
                .aminoAcids()
                .asGroupContainer();
        nonHydrogenAtoms = Selection.on(residueContainer)
                .nonHydrogenAtoms()
                .asAtomContainer();
        nonHydrogenAtoms.atoms()
                .parallel()
                .forEach(this::assignRadius);

        // initialising the sphere points to sample
        spherePoints = generateSpherePoints(numberOfSpherePoints);
        cons = 4.0 * Math.PI / numberOfSpherePoints;

        residueContainer.groups()
                .parallel()
                .forEach(this::assignSingleAsa);
    }

    /**
     * Returns list of 3D coordinates of points on a sphere using the
     * Golden Section Spiral algorithm.
     * @param nSpherePoints the number of points to be used in generating the spherical dot-density
     * @return the collection of generated sphere points
     */
    private List<double[]> generateSpherePoints(int nSpherePoints) {
        List<double[]> points = new ArrayList<>();
        double inc = Math.PI * (3.0 - Math.sqrt(5.0));
        double offset = 2.0 / nSpherePoints;
        for (int k = 0 ; k < nSpherePoints; k++) {
            double y = k * offset - 1.0 + (offset / 2.0);
            double r = Math.sqrt(1.0 - y * y);
            double phi = k * inc;
            points.add(new double[] { Math.cos(phi) * r, y, Math.sin(phi) * r });
        }
        return points;
    }

    /**
     * Gets the van der Waals radius of the given atom following the values defined by Chothia (1976)
     * J.Mol.Biol.105,1-14. NOTE: the vdw values defined by the paper assume no Hydrogens and thus "inflates" slightly
     * the heavy atoms to account for Hydrogens. Thus this method cannot be used in a structure that contains Hydrogens!
     * @param atom the atom whose radius shall be determined
     * @return the atom's radiues
     */
    private double determineRadius(Atom atom) {
        Group parentResidue = atom.getParentGroup();

        switch(atom.getElement()) {
            case H: case D:
                return Element.H.getVDWRadius();
            case O:
                return OXIGEN_VDW;
            case S:
                return SULFUR_VDW;
            case N:
                return atom.getName().equals("NZ") ? TETRAHEDRAL_NITROGEN_VDW : TRIGONAL_NITROGEN_VDW;
            case C:
                String atomName = atom.getName();
                if(atomName.equals("C") || atomName.equals("CE1") || atomName.equals("CE2") || atomName.equals("CE3") ||
                        atomName.equals("CH2") || atomName.equals("CZ") || atomName.equals("CZ2") ||
                        atomName.equals("CZ3")) {
                    return TRIGONAL_CARBON_VDW;
                }
                if (atomName.equals("CA") || atomName.equals("CB") || atomName.equals("CE") || atomName.equals("CG1") ||
                        atomName.equals("CG2")) {
                    return TETRAHEDRAL_CARBON_VDW;
                }
                switch(parentResidue.getGroupInformation().getAminoAcidFamily()) {
                    case PHENYLALANINE: case TRYPTOPHAN: case TYROSINE: case HISTIDINE: case ASPARTIC_ACID: case ASPARAGINE:
                        return TRIGONAL_CARBON_VDW;
                    case PROLINE: case LYSINE: case ARGININE: case METHIONINE: case ISOLEUCINE: case LEUCINE:
                        return TETRAHEDRAL_CARBON_VDW;
                    case GLUTAMIC_ACID: case GLUTAMINE:
                        return atomName.equals("CD") ? TRIGONAL_CARBON_VDW : TETRAHEDRAL_CARBON_VDW;
                    default:
                        //TODO this is not optimal
                        logger.warn("unknown radius for atom of residue  {}", atom.getParentGroup());
                        //TODO what value to return?
                        return TETRAHEDRAL_CARBON_VDW;
                }
            default:
                throw new IllegalArgumentException("unknown case for atom: " + atom);
        }
    }

    /**
     * Returns list of atoms within probe distance to a given atom.
     * @param atom the atom whose neighbor shall be assessed
     * @return all neighbored atoms
     */
    private List<Atom> findNeighbors(Atom atom) {
        final double cutoff = probeSize + probeSize + atom.getFeatureAsDouble(FeatureNames.ATOM_RADIUS);
        return Selection.on(nonHydrogenAtoms)
                .distance(atom, cutoff)
                .asFilteredAtoms()
                .collect(Collectors.toList());
    }

    /**
     * Calculates the accessible surface area (ASA) of an individual atom.
     * @param atom the atom to process
     * @return this atom's ASA
     */
    private double calcSingleAsa(Atom atom) {
        List<Atom> neighborAtoms = findNeighbors(atom);
        double radius = probeSize + atom.getFeatureAsDouble(FeatureNames.ATOM_RADIUS);
        int accessiblePoints = 0;

        for (double[] point : spherePoints) {
            boolean isAccessible = true;
            double[] testPoint = LinearAlgebra3D.add(LinearAlgebra3D.multiply(point, radius), atom.getCoordinates());
            for(Atom neighborAtom : neighborAtoms) {
                double neighborAtomRadius = neighborAtom.getFeatureAsDouble(FeatureNames.ATOM_RADIUS) + probeSize;
                double differenceSquared = LinearAlgebra3D.distanceFast(testPoint, neighborAtom.getCoordinates());
                if (differenceSquared < neighborAtomRadius * neighborAtomRadius) {
                    isAccessible = false;
                    break;
                }
            }
            if (isAccessible) {
                accessiblePoints++;
            }
        }

        return cons * accessiblePoints * radius * radius;
    }
}