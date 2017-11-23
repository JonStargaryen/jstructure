package de.bioforscher.equant;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.net.URL;
import java.nio.channels.Channels;
import java.nio.channels.ReadableByteChannel;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.stream.Stream;

public class EquantConstants {
    private static final Logger logger = LoggerFactory.getLogger(EquantConstants.class);
    public static final Path GIT_DIRECTORY = Paths.get(System.getProperty("user.home")).resolve("git");
    public static final Path PHD_DIRECTORY = GIT_DIRECTORY.resolve("phd_sb_repo");
    public static final Path DATA_DIRECTORY = PHD_DIRECTORY.resolve("data");
    public static final Path EQUANT_DATA_DIRECTORY = DATA_DIRECTORY.resolve("equant");
    public static final Path CASP_DATA_DIRECTORY = EQUANT_DATA_DIRECTORY.resolve("casp");
    public static final Path CASP9_DIRECTORY = CASP_DATA_DIRECTORY.resolve("CASP9");

    public static Stream<String> lines(Path path) {
        try {
            return Files.lines(path)
                    .filter(line -> !line.startsWith("#"));
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    public static Stream<String> lines(String url) {
        try {
            return new BufferedReader(new InputStreamReader(new URL(url).openStream()))
                    .lines()
                    .filter(line -> !line.startsWith("#"));
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    public static Stream<Path> list(Path path) {
        try {
            return Files.list(path);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    public static void write(Path path, String content) {
        ensureDirectoriesExist(path);
        try {
            Files.write(path, content.getBytes());
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    private static void ensureDirectoriesExist(Path path) {
        if(Files.isDirectory(path)) {
            throw new IllegalArgumentException(path.toFile().getAbsolutePath() + " is no regular file - cannot process");
        }

        Path parentDirectory = path.getParent();
        if(!Files.exists(parentDirectory)) {
            logger.info("creating directory '{}'", parentDirectory.toFile().getAbsolutePath());
            try {
                Files.createDirectories(parentDirectory);
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
        }
    }

    public static void move(Path inPath, Path outPath) {
        ensureDirectoriesExist(outPath);
        try {
            Files.move(inPath, outPath);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    public static void move(URL inUrl, Path outPath) {
        ensureDirectoriesExist(outPath);
        try {
            ReadableByteChannel rbc = Channels.newChannel(inUrl.openStream());
            FileOutputStream fos = new FileOutputStream(outPath.toFile());
            fos.getChannel().transferFrom(rbc, 0, Long.MAX_VALUE);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    public static void extractTarGz(Path archivePath) {
        try {
            ProcessBuilder processBuilder = new ProcessBuilder();
            processBuilder.directory(archivePath.getParent().toFile());
            processBuilder.command("tar", "-xzvf", archivePath.toFile().getName());
            processBuilder.start().waitFor();
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        } catch (InterruptedException e) {
            throw new RuntimeException(e);
        }
    }

    public static void delete(Path path) {
        try {
            Files.delete(path);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    public static InputStream newInputStream(Path path) {
        try {
            return Files.newInputStream(path);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }
}
