package design.parser.opm;

import de.bioforscher.jstructure.model.structure.AminoAcid;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Residue;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Describes a trans-membrane helix in OPM terms.
 * Created by S on 29.10.2016.
 */
public class TMHelix {
    private final double tilt;
    private final Residue startResidue;
    private final Residue endResidue;
    private final List<Residue> residues;

    public TMHelix(double tilt, Residue startResidue, Residue endResidue) {
        this.tilt = tilt;
        this.startResidue = startResidue;
        this.endResidue = endResidue;

        int startResNum = startResidue.getResidueNumber();
        int endResNum = endResidue.getResidueNumber();
        Chain chain = startResidue.getParentChain();

        // extract residues
        this.residues = chain.residues()
            .filter(residue -> residue.getResidueNumber() >= startResNum && residue.getResidueNumber() <= endResNum)
            .collect(Collectors.toList());
    }

    public double getTilt() {
        return tilt;
    }

    public Residue getStartResidue() {
        return startResidue;
    }

    public Residue getEndResidue() {
        return endResidue;
    }

    @Override
    public String toString() {
        return getClass().getSimpleName() + " tilt='" + tilt + "°' startResidue='" + startResidue.getPdbName() + "-" +
                startResidue.getResidueNumber() + "' endResidue='" + endResidue.getPdbName() + "-" +
                endResidue.getResidueNumber() + "'";
    }

    public Stream<Residue> residues() {
        return residues.stream();
    }

    public String getSequence() {
        return residues.stream()
                       .map(Residue::getAminoAcid)
                       .map(AminoAcid::getOneLetterCode)
                       .collect(Collectors.joining());
    }
}
