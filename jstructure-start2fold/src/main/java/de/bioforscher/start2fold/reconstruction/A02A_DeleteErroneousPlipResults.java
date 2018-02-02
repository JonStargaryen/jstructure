package de.bioforscher.start2fold.reconstruction;

import de.bioforscher.start2fold.Start2FoldConstants;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.FileVisitOption;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Comparator;

/**
 * Plip contact maps used erroneous constraints - classic.
 */
@Deprecated
public class A02A_DeleteErroneousPlipResults {
    public static void main(String[] args) {
//        deleteContactMaps();

//        deleteReconstructions();

//        deleteTMResults();
    }

    private static void deleteTMResults() {
        Start2FoldConstants.list(Start2FoldConstants.START2FOLD_DIRECTORY.resolve("tmalign"))
                .filter(path -> path.toFile().getName().contains("-plip"))
                .forEach(path -> Start2FoldConstants.list(path)
                        .map(Path::toFile)
                        .forEach(File::delete));
    }

    private static void deleteReconstructions() {
        Start2FoldConstants.list(Start2FoldConstants.START2FOLD_DIRECTORY.resolve("reconstructions"))
                .filter(path -> path.toFile().getName().contains("-plip"))
                .forEach(A02A_DeleteErroneousPlipResults::deleteDirectoryRecursively);
    }

    private static void deleteContactMaps() {
        Start2FoldConstants.list(Start2FoldConstants.START2FOLD_DIRECTORY.resolve("maps"))
                .map(Path::toFile)
                .filter(file -> file.getName().contains("-plip"))
                .forEach(File::delete);
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