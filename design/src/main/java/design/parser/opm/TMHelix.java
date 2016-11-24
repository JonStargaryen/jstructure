package design.parser.opm;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.StructureCollectors;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.selection.Selection;

/**
 * Describes a trans-membrane helix in OPM terms.
 * Created by S on 29.10.2016.
 */
public class TMHelix {
    private final double tilt;
    private final Group startGroup;
    private final Group endGroup;
    private final GroupContainer residues;

    public TMHelix(double tilt, Group startGroup, Group endGroup) {
        this.tilt = tilt;
        this.startGroup = startGroup;
        this.endGroup = endGroup;

        int startResNum = startGroup.getResidueNumber();
        int endResNum = endGroup.getResidueNumber();
        Chain chain = startGroup.getParentChain();

        // extract residues
        this.residues = Selection.on(chain)
                .aminoAcids()
                .asFilteredGroups()
                .filter(residue -> residue.getResidueNumber() >= startResNum && residue.getResidueNumber() <= endResNum)
                .collect(StructureCollectors.toGroupContainer());
    }

    public double getTilt() {
        return tilt;
    }

    public Group getStartGroup() {
        return startGroup;
    }

    public Group getEndGroup() {
        return endGroup;
    }

    @Override
    public String toString() {
        return getClass().getSimpleName() + " tilt='" + tilt + "°' startGroup='" + startGroup.getPdbName() + "-" +
                startGroup.getResidueNumber() + "' endGroup='" + endGroup.getPdbName() + "-" +
                endGroup.getResidueNumber() + "'";
    }

    public GroupContainer getGroupContainer() {
        return residues;
    }

    public String getSequence() {
        return residues.getAminoAcidSequence();
    }
}
