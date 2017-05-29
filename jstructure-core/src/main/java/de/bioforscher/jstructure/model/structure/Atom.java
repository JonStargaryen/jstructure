package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.mathematics.LinearAlgebra;
import de.bioforscher.jstructure.model.Calculable;
import de.bioforscher.jstructure.model.feature.AbstractFeatureable;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.container.StructureContainer;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.Arrays;

/**
 * The most fine-grained element describing a {@link Protein}.
 * Created by S on 27.09.2016.
 */
public class Atom extends AbstractFeatureable implements AtomRecordWriter, CoordinateProvider, StructureContainer,
        Calculable<LinearAlgebra.PrimitiveDoubleArrayLinearAlgebra> {
    private static final Logger logger = LoggerFactory.getLogger(Atom.class);
    public static final float DEFAULT_BFACTOR = 100.0f;
    public static final float DEFAULT_OCCUPANCY = 1.0f;
    public static final String ATOM_PREFIX = "ATOM  ";
    public static final String HETATM_PREFIX = "HETATM";

    private Element element;
    private String name;
    private int pdbSerial;
    private double[] coordinates;
    private Group parentGroup;
    private float occupancy;
    private float bfactor;
    private boolean virtual;
    private String alternativeLocation;
    private String identifier;

    /**
     * Copy constructor.
     * @param atom the original instance
     */
    public Atom(Atom atom) {
        this.element = atom.element;
        this.name = atom.name;
        this.pdbSerial = atom.pdbSerial;
        this.coordinates = Arrays.copyOf(atom.coordinates, 3);
        this.parentGroup = atom.parentGroup;
        this.occupancy = atom.occupancy;
        this.bfactor = atom.bfactor;
        this.virtual = atom.virtual;
        this.alternativeLocation = atom.alternativeLocation;
        this.identifier = atom.identifier;
    }

    //TODO clean-up constructor mess

    /**
     * Constructs an atom.
     * @param name the unique {@link AtomContainer}-wide identifier
     * @param pdbSerial the unique protein-wide identifier
     * @param element this atom's chemical element as Object
     * @param coordinates the spatial coordinates of this atom
     * @param occupancy the occupancy value
     * @param bfactor the b-factor
     * @param alternativeLocation the indicator for alt locations (e.g. A or B) - maybe empty
     */
    public Atom(String name, int pdbSerial, Element element, double[] coordinates, float occupancy, float bfactor, String alternativeLocation) {
        this.name = name;
        this.pdbSerial = pdbSerial;
        this.element = element;
        this.coordinates = coordinates;
        this.occupancy = occupancy;
        this.bfactor = bfactor;
        this.alternativeLocation = alternativeLocation;
    }

    /**
     * Constructs an atom with all information assigned (mostly by being parsed from a <code>PDB</code> file).
     * @param name the unique {@link AtomContainer}-wide identifier
     * @param pdbSerial the unique protein-wide identifier
     * @param elementName this atom's chemical element as <code>String</code> which will be parsed by
     * {@link Element#valueOfIgnoreCase(String)}
     * @param coordinates the spatial coordinates of this atom
     * @param occupancy the occupancy value
     * @param bfactor the b-factor
     */
    public Atom(String name, int pdbSerial, String elementName, double[] coordinates, float occupancy, float bfactor, String alternativeLocation) {
        this(name, pdbSerial, Element.valueOfIgnoreCase(elementName), coordinates, occupancy, bfactor, alternativeLocation);
    }

    /**
     * Constructs an atom without occupancy and bfactor explicitly set (mostly by reconstructing this atom)
     * @param name the unique {@link AtomContainer}-wide identifier
     * @param pdbSerial the u nique protein-wide identifier
     * @param element this atom's chemical element as <code>String</code> which will be parsed by
     * {@link Element#valueOfIgnoreCase(String)}
     * @param coordinates the spatial coordinates of this atom
     */
    public Atom(String name, int pdbSerial, Element element, double[] coordinates) {
        this(name, pdbSerial, element, coordinates, DEFAULT_OCCUPANCY, DEFAULT_BFACTOR, "");
    }

    /**
     * The constructor for a virtual atom with some extra information if necessary.  Calling {@link Atom#isVirtual()}
     * will return <code>true</code>.
     * @param name the name of this atom
     * @param element the element of this atom
     * @param coordinates the coordinates
     */
    public Atom(String name, Element element, double[] coordinates) {
        this(name, 0, element, coordinates, DEFAULT_OCCUPANCY, DEFAULT_BFACTOR, "");
        this.virtual = true;
    }

    /**
     * The constructor for a virtual atom. Calling {@link Atom#isVirtual()} will return <code>true</code>.
     * @param coordinates the coordinates
     */
    public Atom(double[] coordinates) {
        this.coordinates = coordinates;
        this.virtual = true;
    }

    /**
     * Returns a 3D vector of the atom's spatial coordinates.
     * @return a 3D double[]
     */
    @Override
    public double[] getCoordinates() {
        return coordinates;
    }

    /**
     * Returns the <tt>PDB</tt> serial of this atom. A unique way to retrieve atoms.
     * @return an integer
     */
    public int getPdbSerial() {
        return pdbSerial;
    }

    /**
     * A unique identifier for each atom in a {@link AtomContainer}.
     * @return a String, e.g. 'CA' or 'CZ'
     * @see de.bioforscher.jstructure.model.structure.family.AminoAcidFamily#allAtomNames()
     * @see de.bioforscher.jstructure.model.structure.family.AminoAcidFamily#sideChainAtomNames()
     */
    public String getName() {
        return name;
    }


    /**
     * States the type of this atom. Provides indirect access to defined properties of the atom.
     * @return an entry of the {@link Element} enum
     * @see Element
     */
    public Element getElement() {
        return element;
    }

    /**
     * Package-private method to set the parent reference.
     * @param parentGroup the parent
     */
    void setParentGroup(Group parentGroup) {
        this.parentGroup = parentGroup;
    }

    /**
     * Returns the {@link Group} this atom is associated to.
     * @return the parent container
     */
    public Group getParentGroup() {
        return parentGroup != null ? parentGroup : Group.UNKNOWN_GROUP;
    }

    @Override
    public String getPdbRepresentation() {
        try {
            String pdbRecord = AtomRecordProvider.toPDBString(this);
            if (pdbRecord.length() == 0) {
                logger.warn("malformed ATOM record {}", toString());
            }
            return pdbRecord;
        } catch (NullPointerException e) {
            //TODO replace with actual fallback - will fail for virtual atoms without parent references
            return toString();
        }
    }

    @Override
    public String toString() {
        return getClass().getSimpleName() + " '" + getIdentifier() + "' coords='" + Arrays.toString(coordinates);
    }

    @Override
    public String getIdentifier() {
        return identifier == null ? name + "-" + pdbSerial : identifier;
    }

    @Override
    public void setIdentifier(String identifier) {
        this.identifier = identifier;
    }

    /**
     * Returns the occupancy of this atom.
     * @return a float
     */
    public float getOccupancy() {
        return occupancy;
    }

    /**
     * Returns the b-factor of this atom.
     * @return a float
     */
    public float getBfactor() {
        return bfactor;
    }

    /**
     * The indicator for alternative locations - maybe empty if there is none.
     * @return the alt location indicator (if any)
     */
    public String getAlternativeLocation() {
        return alternativeLocation;
    }

    public boolean hasAlternativeLocations() {
        return !alternativeLocation.isEmpty();
    }

    /**
     * Assign new coordinates to this atom.
     * @param coordinates a 3D vector with the coordinates to assign
     */
    @Override
    public void setCoordinates(double[] coordinates) {
        this.coordinates = coordinates;
    }

    public void setBfactor(float bfactor) {
        this.bfactor = bfactor;
    }

    public boolean isVirtual() {
        return virtual;
    }

    public void setPdbSerial(int pdbSerial) {
        this.pdbSerial = pdbSerial;
    }

    @Override
    public LinearAlgebra.PrimitiveDoubleArrayLinearAlgebra calculate() {
        return LinearAlgebra.on(this);
    }
}