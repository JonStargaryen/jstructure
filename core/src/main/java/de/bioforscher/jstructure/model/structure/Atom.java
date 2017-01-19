package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.model.feature.AbstractFeatureContainer;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Locale;

/**
 * The most fine-grained element describing a {@link Protein}.
 * Created by S on 27.09.2016.
 */
public class Atom extends AbstractFeatureContainer implements AtomRecordWriter {
    private final Logger logger = LoggerFactory.getLogger(Atom.class);
    private static final float DEFAULT_BFACTOR = 1.0f;
    private static final float DEFAULT_OCCUPANCY = 1.0f;

    private Element element;
    private String name;
    private int pdbSerial;
    private double[] coordinates;
    private Group parentGroup;
    private float occupancy;
    private float bfactor;
    private boolean virtual;
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
    }

    /**
     * Constructs an atom.
     * @param name the unique {@link AtomContainer}-wide identifier
     * @param pdbSerial the unique protein-wide identifier
     * @param element this atom's chemical element as Object
     * @param coordinates the spatial coordinates of this atom
     * @param occupancy the occupancy value
     * @param bfactor the b-factor
     */
    public Atom(String name, int pdbSerial, Element element, double[] coordinates, float occupancy, float bfactor) {
        this.name = name;
        this.pdbSerial = pdbSerial;
        this.element = element;
        this.coordinates = coordinates;
        this.occupancy = occupancy;
        this.bfactor = bfactor;
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
    public Atom(String name, int pdbSerial, String elementName, double[] coordinates, float occupancy, float bfactor) {
        this(name, pdbSerial, Element.valueOfIgnoreCase(elementName), coordinates, occupancy, bfactor);
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
        this(name, pdbSerial, element, coordinates, DEFAULT_OCCUPANCY, DEFAULT_BFACTOR);
    }

    /**
     * The constructor for a virtual atom with some extra information if necessary.  Calling {@link Atom#isVirtual()}
     * will return <code>true</code>.
     * @param name the name of this atom
     * @param element the element of this atom
     * @param coordinates the coordinates
     */
    public Atom(String name, Element element, double[] coordinates) {
        this(name, 0, element, coordinates, DEFAULT_OCCUPANCY, DEFAULT_BFACTOR);
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

    @Override
    public void setIdentifier(String identifier) {
        this.identifier = identifier;
    }

    /**
     * Returns a 3D vector of the atom's spatial coordinates.
     * @return a 3D double[]
     * @see de.bioforscher.jstructure.mathematics.LinearAlgebra3D
     */
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
    public String composePDBRecord() {
        try {
            String pdbRecord = AtomRecordProvider.toPDBString(this);
            if (pdbRecord.length() == 0) {
                logger.warn("peculiar ATOM record {}", toString());
            }
            return pdbRecord;
        } catch (NullPointerException e) {
            //TODO replace with actual fallback - will fail for virtual atoms without parent references
            return toString();
        }
    }

    @Override
    public String toString() {
        return getClass().getSimpleName() + " identifier='" + getIdentifier() + "' coords='" +
                Arrays.toString(coordinates);
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
     * Assign new coordinates to this atom.
     * @param coordinates a 3D vector with the coordinates to assign
     */
    public void setCoordinates(double[] coordinates) {
        this.coordinates = coordinates;
    }

    public boolean isVirtual() {
        return virtual;
    }

    @Override
    public String getIdentifier() {
        return identifier == null ? name + "-" + pdbSerial : identifier;
    }

    /**
     * Inner class to write <tt>ATOM</tt> records. BioJava-Code.
     */
    static class AtomRecordProvider {
        static DecimalFormat d3 = (DecimalFormat) NumberFormat.getInstance(Locale.US);

        static {
            d3.setMaximumIntegerDigits(4);
            d3.setMinimumFractionDigits(3);
            d3.setMaximumFractionDigits(3);
        }

        static DecimalFormat d2 = (DecimalFormat) NumberFormat.getInstance(Locale.US);

        static {
            d2.setMaximumIntegerDigits(3);
            d2.setMinimumFractionDigits(2);
            d2.setMaximumFractionDigits(2);
        }

        static String toPDBString(Atom atom) {
            Group parentGroup = atom.parentGroup;
            Chain parentChain = parentGroup.getParentChain();
            //TODO check - needed?
            String chainId = parentChain.getChainId() != null ? parentChain.getChainId() : Chain.UNKNOWN_CHAIN.getChainId();
            String record = parentGroup.isLigand() ? "HETATM" : "ATOM  ";
            // format output ...
            String resName = parentGroup.getThreeLetterCode();
            String pdbcode = String.valueOf(parentGroup.getResidueNumber());

            int seri = atom.pdbSerial;
            String serial = String.format("%5d", seri);
            String fullName = formatAtomName(atom);

            Character altLoc = ' ';
            String resseq = String.format("%4s", pdbcode) + " ";

            String x = String.format("%8s", d3.format(atom.coordinates[0]));
            String y = String.format("%8s", d3.format(atom.coordinates[1]));
            String z = String.format("%8s", d3.format(atom.coordinates[2]));
            String occupancy = String.format("%6s", d2.format(atom.occupancy));
            String bfactor = String.format("%6s", d2.format(atom.bfactor));
            String leftResName = String.format("%3s", resName);

            String s = record +
                    serial +
                    " " +
                    fullName +
                    altLoc +
                    leftResName +
                    " " +
                    chainId +
                    resseq +
                    "   " +
                    x +
                    y +
                    z +
                    occupancy +
                    bfactor;

            String e = atom.element.toString().toUpperCase();

            return String.format("%-76s%2s", s, e) + "  ";
        }

        static String formatAtomName(Atom atom) {
            String fullName = null;
            String name = atom.name;
            String element = atom.element.toString().toUpperCase();

            // RULES FOR ATOM NAME PADDING: 4 columns in total: 13, 14, 15, 16

            // if length 4: nothing to do
            if (name.length() == 4) {
                fullName = name;
            } else if (name.length() == 3) {
                fullName = " " + name;
            } else if (name.length() == 2) {
                if (element.equals("C") || element.equals("N") || element.equals("O") || element.equals("P") || element.equals("S")) {
                    fullName = " " + name + " ";
                } else {
                    fullName = name + "  ";
                }
            } else if (name.length() == 1) {
                fullName = " " + name + "  ";
            }

            return fullName;
        }
    }
}