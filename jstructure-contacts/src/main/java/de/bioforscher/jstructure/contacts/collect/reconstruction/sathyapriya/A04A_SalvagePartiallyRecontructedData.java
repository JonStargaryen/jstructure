package de.bioforscher.jstructure.contacts.collect.reconstruction.sathyapriya;

import de.bioforscher.jstructure.contacts.ContactsConstants;
import de.bioforscher.jstructure.contacts.collect.reconstruction.ReconstructionContactMapCreator;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.FileVisitOption;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Comparator;
import java.util.stream.Collectors;

@Deprecated
public class A04A_SalvagePartiallyRecontructedData {
    private static final Path BASE_DIRECTORY = ContactsConstants.RECONSTRUCTION_DIRECTORY.resolve("reconstructions");

    public static void main(String[] args) {
//        deleteOldPlipReconstructions();

        deletePartialResults();

        createMissingConfoldScripts();
    }

    private static void createMissingConfoldScripts() {
        for(ReconstructionContactMapCreator.Sampling sampling : ReconstructionContactMapCreator.Sampling.values()) {
            String samplingPercentage = sampling.name().substring(1);

            String output = createMissingConfoldScripts(sampling.name());
            ContactsConstants.write(BASE_DIRECTORY.getParent().resolve("sampled-" + samplingPercentage + "-reconstruction.sh"),
                    output);
        }
    }

    private static String createMissingConfoldScripts(String samplingName) {
        Path samplingDirectory = BASE_DIRECTORY.resolve(samplingName);
        return ContactsConstants.list(samplingDirectory)
                // filter for directories with missing models
                .filter(path -> {
                    Path stage1 = path.resolve("stage1");
                    return !Files.exists(stage1) || ContactsConstants.list(stage1).noneMatch(model -> model.toFile().getName().contains("_model"));
                })
                .map(path -> {
                    String jobName = path.toFile().getName();
                    String pdbId = jobName.split("-")[0];
                    String mode = jobName.split("-")[1];
                    String mapId = jobName.split("-")[2];
                    return "/home/bittrich/programs/confold_v1.0/confold.pl " +
                            "-seq /home/bittrich/git/phd_sb_repo/data/reconstruction/maps/" + pdbId + "_A.fasta " +
                            "-rr /home/bittrich/git/phd_sb_repo/data/reconstruction/maps/" + samplingName + "/" + pdbId + "_A-" + mode + "-" + mapId + ".rr " +
                            "-o /home/bittrich/git/phd_sb_repo/data/reconstruction/reconstructions/" + samplingName + "/" + jobName + "/";
                })
                .collect(Collectors.joining(System.lineSeparator()));
    }

    private static void deletePartialResults() {
        ContactsConstants.list(BASE_DIRECTORY)
                .filter(path -> !path.toFile().getName().equals("p100"))
                .flatMap(ContactsConstants::list)
                .filter(path -> Files.exists(path.resolve("stage1")))
                .filter(path -> ContactsConstants.list(path.resolve("stage1"))
                        .noneMatch(stage -> stage.toFile().getName().contains("_model")))
                .forEach(A04A_SalvagePartiallyRecontructedData::deleteDirectoryRecursively);
    }

    private static void deleteOldPlipReconstructions() {
        ContactsConstants.list(BASE_DIRECTORY)
                .flatMap(ContactsConstants::list)
                .filter(path -> path.toFile().getName().contains("-plip"))
                .forEach(A04A_SalvagePartiallyRecontructedData::deleteDirectoryRecursively);
    }

    private static void deleteDirectoryRecursively(Path path) {
        try {
            Path input = path.resolve("input");
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
