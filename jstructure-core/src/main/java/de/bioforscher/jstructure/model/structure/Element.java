package de.bioforscher.jstructure.model.structure;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.HashMap;
import java.util.Map;

/**
 * Enumerates all elements which can be represented by an {@link Atom}. Thank BioJava for all these sweet values!
 * Created by S on 28.09.2016.
 */
public enum Element {
    H(1, 1, 1.10f, 1.008f),
    C(6, 2, 1.55f, 12.011f),
    N(7, 2, 1.40f, 14.007f),
    O(8, 2, 1.35f, 16.000f),
    D(1, 1, 1.10f, 1.008f),
    T(1, 1, 1.10f, 1.008f),
    He(2, 1, 2.20f, 4.003f),
    Li(3, 2, 1.22f, 6.940f),
    Be(4, 2, 0.63f, 9.012f),
    B(5, 2, 1.55f, 10.810f),
    F(9, 2, 1.30f, 18.998f),
    Ne(10, 2, 2.02f, 20.170f),
    Na(11, 3, 2.20f, 22.990f),
    Mg(12, 3, 1.50f, 24.305f),
    Al(13, 3, 1.50f, 26.982f),
    Si(14, 3, 2.20f, 28.086f),
    P(15, 3, 1.88f, 30.974f),
    S(16, 3, 1.81f, 32.060f),
    Cl(17, 3, 1.75f, 35.453f),
    Ar(18, 4, 2.77f, 39.948f),
    K(19, 4, 2.39f, 39.102f),
    Ca(20, 4, 1.95f, 40.080f),
    Sc(21, 4, 1.32f, 44.956f),
    Ti(22, 4, 1.95f, 47.880f),
    V(23, 4, 1.06f, 50.040f),
    Cr(24, 4, 1.13f, 51.996f),
    Mn(25, 4, 1.19f, 54.938f),
    Fe(26, 4, 1.95f, 55.847f),
    Co(27, 4, 1.13f, 58.933f),
    Ni(28, 4, 1.24f, 58.710f),
    Cu(29, 4, 1.15f, 63.546f),
    Zn(30, 4, 1.15f, 65.380f),
    Ga(31, 4, 1.55f, 69.720f),
    Ge(32, 4, 2.72f, 72.590f),
    As(33, 4, 0.83f, 74.922f),
    Se(34, 4, 0.90f, 78.960f),
    Br(35, 4, 1.95f, 79.904f),
    Kr(36, 4, 1.90f, 83.800f),
    Rb(37, 5, 2.65f, 85.467f),
    Sr(38, 5, 2.02f, 87.620f),
    Y(39, 5, 1.61f, 88.806f),
    Zr(40, 5, 1.42f, 91.220f),
    Nb(41, 5, 1.33f, 92.906f),
    Mo(42, 5, 1.75f, 95.940f),
    Tc(43, 5, 1.80f, 98.910f),
    Ru(44, 5, 1.20f, 101.070f),
    Rh(45, 5, 1.22f, 102.906f),
    Pd(46, 5, 1.44f, 106.400f),
    Ag(47, 5, 1.55f, 107.868f),
    Cd(48, 5, 1.75f, 112.400f),
    In(49, 5, 1.46f, 114.820f),
    Sn(50, 5, 1.67f, 118.690f),
    Sb(51, 5, 1.12f, 121.750f),
    Te(52, 5, 1.26f, 127.600f),
    I(53, 5, 2.15f, 126.905f),
    Xe(54, 5, 2.10f, 131.300f),
    Cs(55, 6, 3.01f, 132.905f),
    Ba(56, 6, 2.41f, 137.340f),
    La(57, 6, 1.83f, 138.905f),
    Ce(58, 6, 1.86f, 140.120f),
    Pr(59, 6, 1.62f, 140.908f),
    Nd(60, 6, 1.79f, 144.240f),
    Pm(61, 6, 1.76f, 145.000f),
    Sm(62, 6, 1.74f, 150.400f),
    Eu(63, 6, 1.96f, 151.960f),
    Gd(64, 6, 1.69f, 157.250f),
    Tb(65, 6, 1.66f, 158.925f),
    Dy(66, 6, 1.63f, 162.500f),
    Ho(67, 6, 1.61f, 164.930f),
    Er(68, 6, 1.59f, 167.260f),
    Tm(69, 6, 1.57f, 168.934f),
    Yb(70, 6, 1.54f, 173.040f),
    Lu(71, 6, 1.53f, 174.970f),
    Hf(72, 6, 1.40f, 178.490f),
    Ta(73, 6, 1.22f, 180.850f),
    W(74, 6, 1.26f, 183.850f),
    Re(75, 6, 1.30f, 186.200f),
    Os(76, 6, 1.58f, 190.200f),
    Ir(77, 6, 1.22f,  192.220f),
    Pt(78, 6, 1.55f, 195.090f),
    Au(79, 6, 1.45f, 196.967f),
    Hg(80, 6, 1.55f, 200.59f),
    Tl(81, 6, 1.96f, 204.3833f),
    Pb(82, 6, 2.16f, 207.200f),
    Bi(83, 6, 1.73f, 208.981f),
    Po(84, 6, 1.21f, 209.000f),
    At(85, 6, 1.12f, 10.000f),
    Rn(86, 6, 2.30f, 222.000f),
    Fr(87, 7, 3.24f, 223.000f),
    Ra(88, 7, 2.57f, 226.000f),
    Ac(89, 7, 2.12f, 227.000f),
    Th(90, 7, 1.84f, 232.038f),
    Pa(91, 7, 1.60f, 231.036f),
    U(92, 7, 1.75f, 238.029f),
    Np(93, 7, 1.71f, 237.048f),
    Pu(94, 7, 1.67f, 244.000f),
    Am(95, 7, 1.66f, 243.000f),
    Cm(96, 7, 1.65f, 248.000f),
    Bk(97, 7, 1.64f, 247.000f),
    Cf(98, 7, 1.63f, 251.000f),
    Es(99, 7, 1.62f, 254.000f),
    Fm(100, 7, 1.61f, 257.000f),
    Md(101, 7, 1.60f, 256.000f),
    No(102, 7, 1.59f, 254.000f),
    Lr(103, 7, 1.58f, 257.000f),
    /**
     * R-group to represent generic groups that are sometimes present in MDL .sdf
     * files.
     */
    R(104, 0, 0.0f, 0.000f), // this is an R-group
    X(105, 0, 0.0f, 0.000f); // this is an unknown element

    static Logger logger = LoggerFactory.getLogger(Element.class);

    private final int atomicNumber;
    private final int period;
    private final float VDWRadius; // in Angstroms
    private final float atomicMass;
    private static final Map<String, Element> allElements;

    static {
        allElements = new HashMap<>();
        for (Element e : Element.values()){
            allElements.put(e.toString().toLowerCase(), e);
        }
    }

    Element(int atomicNumber, int period, float VDWRadius, float atomicMass) {
        this.atomicNumber = atomicNumber;
        this.period = period;
        this.VDWRadius = VDWRadius;
        this.atomicMass = atomicMass;
    }

    /**
     * Returns the van der Waals radius of this Element.
     * @return the van der Waals radius of this Element, measured in Angstroms.
     */
    public float getVDWRadius() {
        return VDWRadius;
    }

    /**
     * Returns the atomic mass for this Element.
     * @return the atomic mass for this Element, measured in g/mol.
     */
    public float getAtomicMass() {
        return atomicMass;
    }

    /**
     * Returns the Element that corresponds to the specified element symbol. The case
     * of the element symbol is ignored. Example: FE, fe, Fe represent iron.
     * @param elementSymbol element symbol to specify Element.
     * @return the Element specified by the element symbol.
     */
    public static Element resolveElementSymbol(String elementSymbol) throws IllegalArgumentException {
        Element e = allElements.get(elementSymbol.toLowerCase());
        if (e != null) {
            return e;
        }
        logger.warn("Invalid element symbol: {} - falling back to {}",
                elementSymbol,
                Element.X);
        return Element.X;
    }

    public static Element resolveFullAtomName(String atomName, boolean hetAtm) {
        try {
            return resolveFullAtomNameInternal(atomName, hetAtm);
        } catch (IllegalArgumentException e) {
            throw new IllegalArgumentException("could not resolve element from atom name '" + atomName + "'");
        }
    }

    private static Element resolveFullAtomNameInternal(String atomName, boolean hetAtm) {
        String numberFreeAtomName = atomName.replaceAll("\\d","");

        if(atomName.isEmpty()) {
            throw new IllegalArgumentException("could not resolve atom name");
        }

        // fix for Holmium
        if(atomName.startsWith("HO")) {
            return Element.H;
        }

        // ambiguity between CA representing an alpha carbon or calcium
        if(!hetAtm && atomName.startsWith("C")) {
            return Element.C;
        }
        // ambiguity between ND representing an sidechain nitrogen or Nd
        if(!hetAtm && atomName.startsWith("N")) {
            return Element.N;
        }
        Element element = allElements.get(numberFreeAtomName.toLowerCase());
        if(element == null) {
            return resolveFullAtomNameInternal(atomName.substring(0, atomName.length() - 1), hetAtm);
        }
        return element;
    }

    /**
     * Returns <code>true</code> if this Element is Hydrogen.
     * <p>
     * <strong>Note:</strong> Deuterium ({@link #D}) and Tritium ({@link Element#T}) will return
     * <code>true</code> to this method.
     * </p>
     *
     * @return true iff the Element is Hydrogen.
     */
    public boolean isHydrogen() {
        return this == H || this == D || this == T;
    }

    /**
     * Returns <CODE>true</CODE> is the Element is an not Hydrogen (or an
     * isotope of Hydrogen).
     * <p>
     * This method is the exact opposite of {@link #isHydrogen()}.
     * </p>
     *
     * @return true iff Element is no hydrogen atom
     */
    public boolean isHeavyAtom() {
        return !isHydrogen();
    }
}