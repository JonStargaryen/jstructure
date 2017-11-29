package de.bioforscher.jstructure.contacts.collect.reconstruction.start2fold;

import de.bioforscher.jstructure.contacts.ContactsConstants;

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
    private static final Path BASE_DIRECTORY = ContactsConstants.START2FOLD_DIRECTORY.resolve("reconstructions");

    public static void main(String[] args) {
//        deletePartialResults();

        createMissingConventionalSampledConfoldScripts();
    }

    private static void createMissingConventionalSampledConfoldScripts() {
        String output = ContactsConstants.list(BASE_DIRECTORY)
                .filter(path -> path.toFile().getName().contains("-sampled-conventional-"))
                // filter for directories with missing models
                .filter(path -> {
                    Path stage1 = path.resolve("stage1");
                    return !Files.exists(stage1) || ContactsConstants.list(stage1).noneMatch(model -> model.toFile().getName().contains("_model"));
                })
                .map(path -> {
                    String jobName = path.toFile().getName();
                    String pdbId = jobName.split("-")[0];
                    String mapId = jobName.split("-")[3];
                    return "mkdir /home/bittrich/tmp/" + jobName + "/" + System.lineSeparator() +
                            "/home/bittrich/programs/confold_v1.0/confold.pl -seq /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A.fasta -rr /home/bittrich/git/phd_sb_repo/data/reconstruction-start2fold/maps/" + pdbId + "_A-sampled-conventional-" + mapId + ".rr -o /home/bittrich/tmp/" + jobName + "/";
                })
                .collect(Collectors.joining(System.lineSeparator()));
        ContactsConstants.write(BASE_DIRECTORY.getParent().resolve("conventional-sampled-reconstruction.sh"),
                output);
    }

    private static void deletePartialResults() {
        ContactsConstants.list(BASE_DIRECTORY)
                .filter(path -> path.toFile().getName().contains("-sampled-conventional-"))
                .filter(path -> Files.exists(path.resolve("stage1")))
                .filter(path -> ContactsConstants.list(path.resolve("stage1"))
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
