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

    public static AtomBuilder builder(Element element, double[] coordinates) {
        return new AtomBuilder(element, coordinates);
    }

    public static class AtomBuilder {
        Element element;
        String name;
        int pdbSerial = 0;
        double[] coordinates;
        String alternativeLocation = "";
        float occupancy = DEFAULT_OCCUPANCY;
        float bfactor = DEFAULT_BFACTOR;
        boolean virtual;

        AtomBuilder(Element element, double[] coordinates) {
            this.element = element;
            this.coordinates = coordinates;
        }

        public AtomBuilder name(String name) {
            this.name = name;
            return this;
        }

        public AtomBuilder pdbSerial(int pdbSerial) {
            this.pdbSerial = pdbSerial;
            return this;
        }

        public AtomBuilder alternativeLocation(String alternativeLocation) {
            this.alternativeLocation = alternativeLocation;
            return this;
        }

        public AtomBuilder occupancy(float occupancy) {
            this.occupancy = occupancy;
            return this;
        }

        public AtomBuilder bfactor(float bfactor) {
            this.bfactor = bfactor;
            return this;
        }

        public Atom build() {
            this.virtual = pdbSerial == 0;
            if(!virtual && name == null) {
                throw new IllegalArgumentException("no atom name provided for non-virtual atom with pdbSerial " + pdbSerial);
            }
            return new Atom(this);
        }
    }

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

    Atom(AtomBuilder atomBuilder) {
        this.element = atomBuilder.element;
        this.name = atomBuilder.name;
        this.pdbSerial = atomBuilder.pdbSerial;
        this.coordinates = atomBuilder.coordinates;
        this.occupancy = atomBuilder.occupancy;
        this.bfactor = atomBuilder.bfactor;
        this.alternativeLocation = atomBuilder.alternativeLocation;
        this.virtual = atomBuilder.virtual;
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
            e.printStackTrace();
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