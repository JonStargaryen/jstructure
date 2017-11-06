package de.bioforscher.jstructure.membrane.collect.couplings;

import de.bioforscher.jstructure.membrane.MembraneConstants;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.zip.GZIPInputStream;

public class A01_OutputPdbList {
    public static void main(String[] args) {
        MembraneConstants.lines(MembraneConstants.COUPLING_DIRECTORY.resolve("ids.list"))
                .map(line -> line.split("_")[0])
                .distinct()
                .forEach(id -> {
                    try {
                        System.out.println("moving " + id);
                        Path pdbFile = Paths.get("/var/local/pdb/" + id.substring(1, 3) + "/pdb" + id + ".ent.gz");
                        GZIPInputStream inputStream = new GZIPInputStream(Files.newInputStream(pdbFile));
                        Path outDirectory = MembraneConstants.COUPLING_DIRECTORY.resolve("pdb");
                        Files.copy(inputStream, outDirectory.resolve(id + ".pdb"));
                    } catch (IOException e) {
                        throw new UncheckedIOException(e);
                    }
                });
    }
}
