package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.model.structure.aminoacid.NonStandardAminoAcid;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Locale;

import static de.bioforscher.jstructure.model.structure.Atom.ATOM_PREFIX;
import static de.bioforscher.jstructure.model.structure.Atom.HETATM_PREFIX;

/**
 * Class to write <tt>ATOM</tt> records. BioJava-Code.
 * Created by bittrich on 5/23/17.
 */
class AtomRecordProvider {
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
        Group parentGroup = atom.getParentGroup();
        Chain parentChain = parentGroup.getParentChain();
        //TODO check - needed?
        String chainId = parentChain.getChainId() != null ? parentChain.getChainId().getChainId() : Chain.UNKNOWN_CHAIN.getChainId().getChainId();
        String record = isHetAtm(parentGroup) ? HETATM_PREFIX : ATOM_PREFIX;
        // format output ...
        String resName = parentGroup.getThreeLetterCode();
        String pdbcode = String.valueOf(parentGroup.getResidueIdentifier().getResidueNumber());

        int seri = atom.getPdbSerial();
        String serial = String.format("%5d", seri);
        String fullName = formatAtomName(atom);

        String altLoc = atom.hasAlternativeLocations() ? atom.getAlternativeLocation() : " ";
        String resseq;
        if(parentGroup.getResidueIdentifier().hasInsertionCode()) {
            // substring for safety
            resseq = String.format("%4s", pdbcode) + parentGroup.getResidueIdentifier().getInsertionCode().substring(0, 1);
        } else {
            resseq = String.format("%4s", pdbcode) + " ";
        }

        double[] coordinates = atom.getCoordinates();
        String x = String.format("%8s", d3.format(coordinates[0]));
        String y = String.format("%8s", d3.format(coordinates[1]));
        String z = String.format("%8s", d3.format(coordinates[2]));
        String occupancy = String.format("%6s", d2.format(atom.getOccupancy()));
        String bfactor = String.format("%6s", d2.format(atom.getBfactor()));
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

        String e = atom.getElement().toString().toUpperCase();

        return String.format("%-76s%2s", s, e) + "  ";
    }

    private static boolean isHetAtm(Group parentGroup) {
        //TODO is this check for HETATM annotation correct?
        return parentGroup.isLigand() || parentGroup instanceof NonStandardAminoAcid;
    }

    private static String formatAtomName(Atom atom) {
        String fullName = null;
        String name = atom.getName();
        String element = atom.getElement().toString().toUpperCase();

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
