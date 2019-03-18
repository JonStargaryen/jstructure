package de.bioforscher.aars.motifs;

import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.StringJoiner;

/**
 * For a small number of sequences to computation failed. Generate the sane list which should not result in any
 * exceptions when accessing the results.
 */
public class S03_CreateSaneNrPdbList {
    public static void main(String[] args) throws IOException {
        StringJoiner output = new StringJoiner(System.lineSeparator());
        Path motifPath = Paths.get("/home/bittrich/git/generic-motif-search/data/extracted-motifs/");

        Files.list(motifPath)
                // chains are valid when they contain a matchable motif with more than 1 residue
                //TODO potentially check for spatial proximity
                //TODO potentially check for present backbone atoms
                .filter(S03_CreateSaneNrPdbList::filterValidChains)
                .map(Path::toFile)
                .map(File::getName)
                .map(name -> name.split("\\.")[0])
                .map(name -> name.split("_"))
                .map(split -> split[0] + "\t" + split[1])
                .sorted()
                .forEach(output::add);

        Files.write(Paths.get("/home/bittrich/git/generic-motif-search/data/nrpdb-filtered.list"),
                output.toString().getBytes());
    }

    private static boolean filterValidChains(Path path) {
        System.out.println("processing " + path);
        Structure structure = StructureParser.fromPath(path).parse();
        return structure.getGroups().size() > 1;
    }
}
