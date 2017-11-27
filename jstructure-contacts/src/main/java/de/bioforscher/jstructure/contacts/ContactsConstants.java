package de.bioforscher.jstructure.contacts;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.testutil.FileUtils;

import java.nio.file.Path;

public class ContactsConstants extends FileUtils {
    public static final Path RECONSTRUCTION_DIRECTORY = DATA_DIRECTORY.resolve("reconstruction");
    public static final Path START2FOLD_DIRECTORY = DATA_DIRECTORY.resolve("reconstruction-start2fold");

    public static Atom getBetaCarbon(AminoAcid aminoAcid) {
        return aminoAcid.atoms()
                .filter(atom -> atom.getName().equals("CB"))
                .findFirst()
                .orElse(aminoAcid.getCa());
    }
}
