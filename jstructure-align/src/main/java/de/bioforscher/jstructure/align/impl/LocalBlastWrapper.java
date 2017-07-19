package de.bioforscher.jstructure.align.impl;

import de.bioforscher.jstructure.align.AlignmentException;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Wrapper for a local blast executable.
 * Created by bittrich on 7/18/17.
 */
public class LocalBlastWrapper {
    /**
     * Execute a PSI-BLAST run using the <code>nr</code> database.
     *
     * <pre>original: blastpgp -d [nr database] -i [target FASTA sequence] -e 1e-10 -J t -u 1 -j 3
     * 'translated': psiblast -db nr -query /path/file -evalue 1e-10 -num_iterations 3 -outfmt 6</pre>
     */
    public PsiBlastResult executePsiBlastRun(String fastaSequence) {
        try {
            Path tmpDirectory = Files.createTempDirectory("psiblast");
            Path sequencePath = tmpDirectory.resolve("query.fsa");
            Path outputPath = tmpDirectory.resolve("output.tsv");
            Path matrixPath = tmpDirectory.resolve("pssm.mtx");

            Files.write(sequencePath, fastaSequence.getBytes());

            ProcessBuilder processBuilder = new ProcessBuilder("psiblast",
                    "-db",
                    "/var/local/blastdb/swissprot",
                    "-query",
                    sequencePath.toFile().getAbsolutePath(),
                    "-evalue",
                    "1e-10",
                    "-num_iterations",
                    "3",
                    "-num_threads",
                    String.valueOf(Math.max((int) (0.5 * Runtime.getRuntime().availableProcessors()), 1)),
                    "-outfmt",
                    "6",
                    "-out",
                    outputPath.toFile().getAbsolutePath(),
                    "-out_ascii_pssm",
                    matrixPath.toFile().getAbsolutePath());

            processBuilder.start().waitFor();
            PsiBlastResult result = new PsiBlastResult(parseResultFile(outputPath), parseMatrixFile(matrixPath));

            Files.delete(sequencePath);
            Files.delete(outputPath);
            Files.delete(matrixPath);

            return result;
        } catch (Exception e) {
            throw new AlignmentException(e);
        }
    }

    List<String> parseResultFile(Path file) throws IOException {
        return parseResultFile(Files.lines(file));
    }

    List<String> parseResultFile(Stream<String> fileStream) {
        return fileStream
                .map(line -> line.split("\\s+"))
                .filter(split -> split.length == 12)
                .map(split -> split[1])
                .map(accession -> accession.split("\\.")[0])
                .distinct()
                .collect(Collectors.toList());
    }

    List<double[]> parseMatrixFile(Path file) throws IOException {
        return parseMatrixFile(Files.lines(file));
    }

    List<double[]> parseMatrixFile(Stream<String> fileStream) {
        return fileStream.filter(line -> line.length() == 181)
                .map(line -> line.split("\\s+"))
                .map(split -> Stream.concat(IntStream.range(3, 23).boxed(), Stream.of(43))
                        .mapToDouble(position -> Double.valueOf(split[position]))
                        .toArray())
                .collect(Collectors.toList());
    }

    public static class PsiBlastResult {
        private static final String AMINO_ACID_HEADER = "A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V";
        private static final Pattern PATTERN = Pattern.compile("\\s+");
        private static final List<AminoAcid.Family> KEY_SET = PATTERN.splitAsStream(AMINO_ACID_HEADER)
                .map(AminoAcid.Family::resolveOneLetterCode)
                .collect(Collectors.toList());
        private final List<String> accessions;
        private final List<Double> information;
        private final List<Map<AminoAcid.Family, Double>> exchanges;

        public PsiBlastResult(List<String> accessions, List<double[]> rawPssm) {
            this.accessions = accessions;
            this.information = rawPssm.stream()
                    .map(array -> array[20])
                    .collect(Collectors.toList());
            this.exchanges = new ArrayList<>();
            for(double[] array : rawPssm) {
                Map<AminoAcid.Family, Double> line = new EnumMap<>(AminoAcid.Family.class);
                for(int i = 0; i < 20; i++) {
                    line.put(KEY_SET.get(i), array[i]);
                }
                exchanges.add(line);
            }
        }

        public List<String> getAccessions() {
            return accessions;
        }

        public List<Double> getInformation() {
            return information;
        }

        public List<Map<AminoAcid.Family, Double>> getExchanges() {
            return exchanges;
        }
    }
}
