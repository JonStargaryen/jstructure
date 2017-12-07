package de.bioforscher.jstructure.membrane.collect.foldingcores;

import de.bioforscher.jstructure.membrane.MembraneConstants;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

public class A01_ExtractSequences {
    public static void main(String[] args) {
        MembraneConstants.list(MembraneConstants.FOLDING_CORES_DIRECTORY.resolve("start2fold"))
                .forEach(A01_ExtractSequences::handleFile);
    }

    private static void handleFile(Path path) {
        try {
            String id = Files.lines(path)
                    .peek(System.out::println)
                    .filter(line -> line.startsWith("#pdb:"))
                    .findFirst()
                    .get()
                    .split(": ")[1];
            String pdbId = id.split("_")[0];
            String chainId = id.split("_")[1];

            Structure structure = StructureParser.source(pdbId).parse();
            Chain chain = structure.select()
                    .chainId(chainId)
                    .asChain();

            MembraneConstants.write(MembraneConstants.FOLDING_CORES_DIRECTORY.resolve("couplings").resolve(pdbId + "_" + chainId + ".fasta"),
                    ">" + pdbId + "_" + chainId + System.lineSeparator() +
                            chain.getAminoAcidSequence());
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
