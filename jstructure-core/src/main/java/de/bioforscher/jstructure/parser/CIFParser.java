package de.bioforscher.jstructure.parser;

import de.bioforscher.jstructure.model.structure.family.GroupInformation;

import java.io.*;
import java.net.URL;
import java.util.HashMap;
import java.util.Map;

/**
 * Fetches and parses ligand information from the PDB.
 * TODO one should disable this by a flag as every operation will require an active internet connection
 * Created by S on 05.01.2017.
 */
public class CIFParser {
    /**
     * The URL which can be used to fetch ligand information.
     */
    private static final String CIF_FETCH_URL = "https://files.rcsb.org/ligands/view/%s.cif";
    /**
     * Collection of already parsed group information.
     */
    private static final Map<String, GroupInformation> knownGroups;

    static {
        knownGroups = new HashMap<>();
        knownGroups.put("UNK", GroupInformation.UNKNOWN_AMINO_ACID);
    }

    public static GroupInformation parseLigandInformation(String ligandId) {
        if(knownGroups.containsKey(ligandId)) {
            return knownGroups.get(ligandId);
        }

        try {
            GroupInformation type = parseLigandInformation(new URL(String.format(CIF_FETCH_URL, ligandId)).openStream());
            knownGroups.put(ligandId, type);
            return type;
        } catch (IOException e) {
            throw new UncheckedIOException("no correct ligand identifier or no internet connection available", e);
        }
    }

    private static GroupInformation parseLigandInformation(InputStream inputStream) throws IOException {
        try (InputStreamReader inputStreamReader = new InputStreamReader(inputStream)) {
            try (BufferedReader bufferedReader = new BufferedReader(inputStreamReader)) {
                return GroupInformation.builder().parseLines(bufferedReader.lines());
            }
        }
    }
}
