package de.bioforscher.jstructure.efr.analysis;

import de.bioforscher.testutil.FileUtils;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Compute EF probability of each sequence position in the representative sequences.
 */
public class A01_ExtractAARSSequences {
    private static final Path SEQUENCE_INPUT = Paths.get("/home/bittrich/git/aars_data/T04_representative_sequences/");
    private static final Path SEQUENCE_OUTPUT = Paths.get("/home/bittrich/git/phd_sb_repo/data/aars/sequences/");
    private static final List<String> CLASS_NAMES = Stream.of("C1", "C2")
            .collect(Collectors.toList());

    public static void main(String[] args) throws IOException {
        StringJoiner efoldmineScript = new StringJoiner(System.lineSeparator());
        for(String className : CLASS_NAMES) {
            List<String> lines = Files.readAllLines(SEQUENCE_INPUT.resolve(className + "_representatives_cluster.fasta"));
            String lastId = "";
            for (String line : lines) {
                if (line.startsWith(">")) {
                    lastId = line.substring(1);
                } else {
                    FileUtils.write(SEQUENCE_OUTPUT.resolve(lastId + ".fasta"),
                            ">" + lastId + System.lineSeparator() + line);
                    efoldmineScript.add("python2 EFoldMine.py sequences/" + lastId + ".fasta -o ef-predictions/" + lastId + ".out");
                }
            }
        }

        FileUtils.write(SEQUENCE_OUTPUT.getParent().resolve("efoldmine.sh"), efoldmineScript.toString());
    }
}
