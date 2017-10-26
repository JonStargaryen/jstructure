package de.bioforscher.jstructure.membrane.modularity.division;

import de.bioforscher.jstructure.membrane.MembraneConstants;
import de.bioforscher.jstructure.model.structure.StructureParser;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.stream.Stream;

public class A01_WriteFastaSequences {
    public static void main(String[] args) {
        Stream.of(MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY,
                MembraneConstants.PDBTM_NR_BETA_DATASET_DIRECTORY)
                .map(path -> path.resolve("plip"))
                .forEach(A01_WriteFastaSequences::handleDirectory);
    }

    private static void handleDirectory(Path path) {
        MembraneConstants.list(path)
                .map(Path::toFile)
                .map(File::getName)
                .map(id -> id.split("\\.")[0])
                .forEach(A01_WriteFastaSequences::handleChainId);
    }

    private static void handleChainId(String id) {
        System.out.println(id);
        String pdbId = id.split("_") [0];
        String chainId = id.split("_")[1];

        String sequence = ">" + id + System.lineSeparator() + StructureParser.source(pdbId)
                .minimalParsing(true)
                .parse()
                .select()
                .chainName(chainId)
                .asChain()
                .getAminoAcidSequence();

        MembraneConstants.write(Paths.get("/home/bittrich/Downloads/dynamine/dynamine/predictor/" + id + ".fasta"), sequence);
    }
}
