package de.bioforscher.jstructure.model.structure.family;

/**
 * Represents a group of atomic or molecular structures such as ions, amino acids or nucleotides.
 * Created by S on 05.01.2017.
 */
public interface AtomicFamily {
    /**
     * The one-letter code representing this group.
     * @return the one-letter code
     */
    String getOneLetterCode();

    /**
     * The three-letter code representing this group.
     * @return the (up-to-)three-letter-code
     */
    String getThreeLetterCode();
}
