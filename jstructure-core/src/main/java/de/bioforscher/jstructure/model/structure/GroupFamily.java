package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.model.structure.identifier.ResidueIdentifier;

/**
 * Specifies properties of group families such as {@link de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid.Family}.
 * Created by bittrich on 7/12/17.
 */
public interface GroupFamily {
    Class<? extends Group> getRepresentingClass();

    GroupPrototype getGroupPrototype();

    Group createGroup(String pdbName, ResidueIdentifier residueIdentifier, boolean ligand);

    String getThreeLetterCode();
}
