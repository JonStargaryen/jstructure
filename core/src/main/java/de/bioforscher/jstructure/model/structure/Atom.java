package de.bioforscher.jstructure.model.structure;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;

/**
 * The most fine-grained element describing a {@link Protein}.
 * Created by S on 27.09.2016.
 */
public class Atom implements AtomRecordWriter, Container {
    public static final float DEFAULT_BFACTOR = 1.0f;
    public static final float DEFAULT_OCCUPANCY = 1.0f;
    private Element element;
    private String name;
    private int pdbSerial;
    private double[] coordinates;
    private Group parentGroup;
    private float occupancy;
    private float bfactor;
    private boolean virtual;
    private Map<Enum, Object> featureMap;

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
        this.featureMap = new HashMap<>();
    }

    /**
     * Constructs an atom.
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
     *
     * @param name
     * @param pdbSerial
     * @param element
     * @param coordinates
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
        this(name, Integer.MIN_VALUE, element, coordinates, DEFAULT_OCCUPANCY, DEFAULT_BFACTOR);
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
     * @see AminoAcid#allAtomNames()
     * @see AminoAcid#sideChainAtomNames()
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
     * Returns the {@link Residue} this atom is associated to.
     * @return the parent container
     */
    public Group getParentGroup() {
        return parentGroup;
    }

    @Override
    public String composePDBRecord() {
        String pdbRecord = AtomRecordProvider.toPDBString(this);
        if (pdbRecord.length() == 0) {
            //TODO logging/error-handling
            System.out.println(this);
        }
        return pdbRecord;
    }

    @Override
    public String toString() {
        return this.getClass().getSimpleName() + " name='" + this.name + "' coords='" +
                Arrays.toString(this.coordinates) + "' element='" + this.element + "'";
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

    @Override
    public Map<Enum, Object> getFeatureMap() {
        return featureMap;
    }

    public boolean isVirtual() {
        return virtual;
    }

    /**
     * Inner class to write <tt>ATOM</tt> records. BioJava-Code.
     */
    static final class AtomRecordProvider {
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

        public static String toPDBString(Atom atom) {
            Group parentGroup = atom.parentGroup;
            String record = parentGroup instanceof Residue ? "ATOM  " : "HETATM";
            // format output ...
            String resName = parentGroup.getPdbName();
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
            String tempfactor = String.format("%6s", d2.format(atom.bfactor));
            String leftResName = String.format("%3s", resName);

            StringBuilder s = new StringBuilder();
            s.append(record);
            s.append(serial);
            s.append(" ");
            s.append(fullName);
            s.append(altLoc);
            s.append(leftResName);
            s.append(" ");
            s.append(parentGroup.getParentChain().getChainId());
            s.append(resseq);
            s.append("   ");
            s.append(x);
            s.append(y);
            s.append(z);
            s.append(occupancy);
            s.append(tempfactor);

            String e = atom.element.toString().toUpperCase();

            return String.format("%-76s%2s", s.toString(), e) + "  ";
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