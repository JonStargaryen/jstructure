package de.bioforscher.jstructure.model.structure;

/**
 * Specifies properties of group families such as {@link de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid.Family}.
 * Created by bittrich on 7/12/17.
 */
public interface GroupFamily {
    Class<? extends Group> getRepresentingClass();

    GroupPrototype getGroupPrototype();

    String getThreeLetterCode();
}
