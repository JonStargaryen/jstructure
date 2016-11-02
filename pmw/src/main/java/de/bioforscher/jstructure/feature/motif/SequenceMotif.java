package de.bioforscher.jstructure.feature.motif;

import de.bioforscher.jstructure.model.structure.AminoAcid;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Residue;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * The container object of a found sequence motif.
 * Created by S on 02.10.2016.
 */
public class SequenceMotif {
    private final SequenceMotifDefinition motifDefinition;
    private final Residue startResidue;
    private final Residue endResidue;
    private final List<Residue> residues;

    public SequenceMotif(SequenceMotifDefinition candidate, Residue startResidue, Residue endResidue) {
        this.motifDefinition = candidate;
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

    public Stream<Residue> residues() {
        return residues.stream();
    }

    public String getSequence() {
        return residues.stream()
                .map(Residue::getAminoAcid)
                .map(AminoAcid::getOneLetterCode)
                .collect(Collectors.joining());
    }

    public Residue getEndResidue() {
        return endResidue;
    }

    public Residue getStartResidue() {
        return startResidue;
    }

    public SequenceMotifDefinition getMotifDefinition() {
        return motifDefinition;
    }

    @Override
    public String toString() {
        return getClass().getSimpleName() + " motifDefinition='" + motifDefinition + "' sequence='" + getSequence() +
                "' startResidue='" + startResidue.getPdbName() + "-" + startResidue.getResidueNumber() +
                "' endResidue='" + endResidue.getPdbName() + "-" + endResidue.getResidueNumber() + "'";
    }
}