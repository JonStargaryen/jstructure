package de.bioforscher.jstructure.efr;

import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.testutil.FileUtils;

import java.nio.file.Path;
import java.util.List;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class Start2FoldConstants extends FileUtils {
    public static final String START2FOLD_URL = "http://start2fold.eu/";
    public static final String START2FOLD_ENTRIES_URL = "http://start2fold.eu/ids";
    public static final String START2FOLD_PROTEINS_URL = "http://start2fold.eu/proteins";

    public static final Path BASE_DIRECTORY = FileUtils.DATA_DIRECTORY.resolve("start2fold");
    public static final Path XML_DIRECTORY = BASE_DIRECTORY.resolve("xml");
    public static final Path FASTA_DIRECTORY = BASE_DIRECTORY.resolve("fasta");
    public static final Path PDB_DIRECTORY = BASE_DIRECTORY.resolve("pdb");
    public static final Path COUPLING_DIRECTORY = BASE_DIRECTORY.resolve("coupling");
    public static final Path PANCSA_LIST = BASE_DIRECTORY.resolve("pancsa.list");
    public static final Path PYMOL_DIRECTORY = BASE_DIRECTORY.resolve("pymol");
    public static final Path EQUANT_DIRECTORY = BASE_DIRECTORY.resolve("equant");
    public static final Path STATISTICS_DIRECTORY = BASE_DIRECTORY.resolve("statistics");
    public static final Path START2FOLD_DIRECTORY = DATA_DIRECTORY.resolve("reconstructions-start2fold");
    public static final Path THERMOSTABILITY_DIRECTORY = DATA_DIRECTORY.resolve("thermostability");

    public static List<Integer> extractFunctionalResidueNumbers(String[] split) {
        return Pattern.compile(",")
                .splitAsStream(split[5].replaceAll("\\[", "").replaceAll("]", ""))
                .flatMapToInt(numberString -> {
                    if(!numberString.contains("-")) {
                        return IntStream.of(Integer.valueOf(numberString));
                    }
                    String[] numberStringSplit = numberString.split("-");
                    return IntStream.rangeClosed(Integer.valueOf(numberStringSplit[0]),
                            Integer.valueOf(numberStringSplit[1]));
                })
                .boxed()
                .collect(Collectors.toList());
    }

    public static TMAlignResult parseTMAlignResultFile(Path tmAlignPath) {
        return new TMAlignResult(tmAlignPath);
    }

    public static boolean areNonCovalentGroups(Group group1, Group group2) {
        return Math.abs(group1.getResidueIndex() - group2.getResidueIndex()) > 1;
    }

    public static class TMAlignResult {
        private final double rmsd;
        private final double tmscore;

        TMAlignResult(Path tmAlignPath) {
            List<String> lines;
            try(Stream<String> stream = Start2FoldConstants.lines(tmAlignPath)) {
                lines = stream.collect(Collectors.toList());
            }
            double rmsd = -1;
            double tmscore = -1;
            for(String line : lines) {
                if(line.startsWith("Aligned length=")) {
                    rmsd = Double.valueOf(line.split("RMSD=")[1].trim().split(",")[0].trim());
                }
                if(line.startsWith("TM-score=")) {
                    tmscore = Double.valueOf(line.split("TM-score=")[1].trim().split("\\(")[0].trim());
                }
            }

            if(rmsd == -1 || tmscore == -1) {
                throw new IllegalArgumentException("tmalign result file malformed - " + tmAlignPath.toFile().getAbsolutePath());
            }

            this.rmsd = rmsd;
            this.tmscore = tmscore;
        }

        public double getRmsd() {
            return rmsd;
        }

        public double getTmscore() {
            return tmscore;
        }
    }
}
