package de.bioforscher.jstructure.contacts.collect.reconstruction.sathyapriya;

import de.bioforscher.jstructure.contacts.ContactsConstants;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;

/**
 * Segmentation fault during tmalign call for 1e6k. Cannot be fixed by rewriting the PDB file.
 */
@Deprecated
public class A05A_FixPdbStructure {
    public static void main(String[] args) {
        Structure structure = StructureParser.source("1e6k").parse();
        ContactsConstants.write(ContactsConstants.RECONSTRUCTION_DIRECTORY.resolve("pdb").resolve("1e6k.pdb"),
                structure.getPdbRepresentation());
    }
}
