package de.bioforscher.jstructure.align.impl;

import de.bioforscher.jstructure.align.AlignmentException;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.feature.SingleValueFeatureContainerEntry;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.stream.Collectors;
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
//                    "/var/local/blastdb/nr",
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

    List<Double> parseMatrixFile(Path file) throws IOException {
        return parseMatrixFile(Files.lines(file));
    }

    List<Double> parseMatrixFile(Stream<String> fileStream) {
        return fileStream.filter(line -> line.length() == 181)
                .map(line -> line.split("\\s+"))
                .map(split -> split[43])
                .map(Double::valueOf)
                .collect(Collectors.toList());
    }

    public static class PsiBlastResult {
        private final List<String> accessions;
        private final List<Double> conservation;

        public PsiBlastResult(List<String> accessions, List<Double> conservation) {
            this.accessions = accessions;
            this.conservation = conservation;
        }

        public List<String> getAccessions() {
            return accessions;
        }

        public List<Double> getConservation() {
            return conservation;
        }
    }

    public static class PSSMConservationScore extends FeatureContainerEntry implements SingleValueFeatureContainerEntry<Double> {
        private final double score;

        public PSSMConservationScore(double score) {
            super(null);
            this.score = score;
        }

        @Override
        public Double getValue() {
            return score;
        }

        public double getScore() {
            return score;
        }
    }
}
