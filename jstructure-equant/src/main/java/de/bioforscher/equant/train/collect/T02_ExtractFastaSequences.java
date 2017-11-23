package de.bioforscher.equant.train.collect;

import de.bioforscher.equant.EquantConstants;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;

import java.nio.file.Path;

/**
 * Extract FASTA sequences for further analysis.
 */
public class T02_ExtractFastaSequences {
    private static final Path DIRECTORY = EquantConstants.CASP9_DIRECTORY;
    private static final Path OUTPUT_DIRECTORY = DIRECTORY.resolve("fasta");

    public static void main(String[] args) {
        EquantConstants.list(DIRECTORY.resolve("targets"))
                .forEach(T02_ExtractFastaSequences::extractFastaSequence);
    }

    private static void extractFastaSequence(Path path) {
        String name = path.toFile().getName().split("\\.")[0];
        Structure structure = StructureParser.source(path).parse();
        String sequence = structure.getAminoAcidSequence();
        String header = ">" + name + System.lineSeparator();
        EquantConstants.write(OUTPUT_DIRECTORY.resolve(name + ".fasta"), header + sequence);
    }
}
