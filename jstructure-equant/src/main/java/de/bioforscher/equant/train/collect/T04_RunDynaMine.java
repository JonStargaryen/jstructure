package de.bioforscher.equant.train.collect;

import de.bioforscher.equant.EquantConstants;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Optional;

/**
 * Processes all FASTA sequences of targets and runs DynaMine.
 * Manually set up all FASTA sequences in /equant/dynamine.
 * Execute run.sh.
 * This class will merely clean-up and move the results.
 */
public class T04_RunDynaMine {
    private static final Path DIRECTORY = EquantConstants.CASP9_DIRECTORY;

    public static void main(String[] args) {
        Path dynamineDir = EquantConstants.EQUANT_DATA_DIRECTORY.resolve("dynamine");

        // delete fasta files
        EquantConstants.list(dynamineDir)
                .filter(path -> !Files.isDirectory(path))
                .filter(path -> path.toFile().getName().endsWith(".fasta"))
                .forEach(EquantConstants::delete);

        // move results
        EquantConstants.list(dynamineDir)
                .filter(Files::isDirectory)
                .filter(path -> !path.toFile().getName().equals("predictor"))
                .map(path -> EquantConstants.list(path).findFirst())
                .filter(Optional::isPresent)
                .map(Optional::get)
                .map(path ->  EquantConstants.list(path)
                        .filter(candidate -> candidate.toFile().getName().endsWith(".pred"))
                        .findFirst())
                .filter(Optional::isPresent)
                .map(Optional::get)
                .forEach(path -> EquantConstants.move(path,
                        DIRECTORY.resolve("dynamine")
                                .resolve(path.getParent().getParent().toFile().getName() + ".dynamine")));

        // delete result folders
        EquantConstants.list(dynamineDir)
                .filter(Files::isDirectory)
                .filter(path -> !path.toFile().getName().equals("predictor"))
                .forEach(EquantConstants::delete);
    }
}
