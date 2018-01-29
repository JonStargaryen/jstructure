package de.bioforscher.start2fold.reconstruction;

import de.bioforscher.start2fold.Start2FoldConstants;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.FileVisitOption;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Comparator;
import java.util.stream.Collectors;

/**
 * CNS cannot handle long file paths. Of course... Limit is 132 and easy workaround is to write results to folder with
 * shorter name. Happened only for sampled conventional.
 */
@Deprecated
public class A03A_SalvagePartiallyRecontructedData {
    private static final Path BASE_DIRECTORY = Start2FoldConstants.START2FOLD_DIRECTORY.resolve("reconstructions");

    public static void main(String[] args) {
//        deletePartialResults();

//        createMissingConventionalSampledConfoldScripts();

        createMissingFullAndDeletedConfoldScripts();
    }

    private static void createMissingFullAndDeletedConfoldScripts() {
        String fullOutput = Start2FoldConstants.list(Start2FoldConstants.START2FOLD_DIRECTORY.resolve("maps"))
                .filter(path -> path.toFile().getName().endsWith("-full-conventional.rr"))
                .map(path -> {
                    String jobName = path.toFile().getName().split("\\.")[0];
                    String pdbId = jobName.split("-")[0];
                    String mapId = "1";
                    return "mkdir /home/bittrich/tmp/" + jobName + "/" + System.lineSeparator() +
                            "/home/bittrich/programs/confold_v1.0/confold.pl " +
                            "-seq /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + ".fasta " +
                            "-rr /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "-full-conventional.rr " +
                            "-o /home/bittrich/tmp/" + jobName + "/";
                })
                .collect(Collectors.joining(System.lineSeparator()));
        Start2FoldConstants.write(BASE_DIRECTORY.getParent().resolve("conventional-full-reconstruction.sh"),
                fullOutput);

        String deletedOutput = Start2FoldConstants.list(Start2FoldConstants.START2FOLD_DIRECTORY.resolve("maps"))
                .filter(path -> path.toFile().getName().endsWith("-deleted-conventional.rr"))
                .map(path -> {
                    String jobName = path.toFile().getName().split("\\.")[0];
                    String pdbId = jobName.split("-")[0];
                    String mapId = "1";
                    return "mkdir /home/bittrich/tmp/" + jobName + "/" + System.lineSeparator() +
                            "/home/bittrich/programs/confold_v1.0/confold.pl " +
                            "-seq /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + ".fasta " +
                            "-rr /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "-deleted-conventional.rr " +
                            "-o /home/bittrich/tmp/" + jobName + "/";
                })
                .collect(Collectors.joining(System.lineSeparator()));
        Start2FoldConstants.write(BASE_DIRECTORY.getParent().resolve("conventional-deleted-reconstruction.sh"),
                deletedOutput);
    }

    private static void createMissingConventionalSampledConfoldScripts() {
        String output = Start2FoldConstants.list(BASE_DIRECTORY)
                .filter(path -> path.toFile().getName().contains("-sampled-conventional-"))
                // filter for directories with missing models
                .filter(path -> {
                    Path stage1 = path.resolve("stage1");
                    return !Files.exists(stage1) || Start2FoldConstants.list(stage1).noneMatch(model -> model.toFile().getName().contains("_model"));
                })
                .map(path -> {
                    String jobName = path.toFile().getName();
                    String pdbId = jobName.split("-")[0];
                    String mapId = jobName.split("-")[3];
                    return "mkdir /home/bittrich/tmp/" + jobName + "/" + System.lineSeparator() +
                            "/home/bittrich/programs/confold_v1.0/confold.pl -seq /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A.fasta -rr /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A-sampled-conventional-" + mapId + ".rr -o /home/bittrich/tmp/" + jobName + "/";
                })
                .collect(Collectors.joining(System.lineSeparator()));
        Start2FoldConstants.write(BASE_DIRECTORY.getParent().resolve("conventional-sampled-reconstruction.sh"),
                output);
    }

    private static void deletePartialResults() {
        Start2FoldConstants.list(BASE_DIRECTORY)
                .filter(path -> path.toFile().getName().contains("-sampled-conventional-"))
                .filter(path -> Files.exists(path.resolve("stage1")))
                .filter(path -> Start2FoldConstants.list(path.resolve("stage1"))
                        .noneMatch(stage -> stage.toFile().getName().contains("_model")))
                .forEach(A03A_SalvagePartiallyRecontructedData::deleteDirectoryRecursively);
    }

    private static void deleteDirectoryRecursively(Path path) {
        try {
            Path input = path.resolve("input");
            System.out.println("about to delete: " + input);
            if (Files.exists(input)) {
                Files.walk(input, FileVisitOption.FOLLOW_LINKS)
                        .sorted(Comparator.reverseOrder())
                        .map(Path::toFile)
                        .peek(System.out::println)
                        .forEach(File::delete);
            }
            Path stage1 = path.resolve("stage1");
            if (Files.exists(stage1)) {
                Files.walk(stage1, FileVisitOption.FOLLOW_LINKS)
                        .sorted(Comparator.reverseOrder())
                        .map(Path::toFile)
                        .peek(System.out::println)
                        .forEach(File::delete);
            }
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }
}
