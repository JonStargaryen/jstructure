package bioforscher.aggregation;

import bioforscher.Constants;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.net.URL;
import java.nio.channels.Channels;
import java.nio.channels.ReadableByteChannel;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

/**
 * Run kink finder on data set.
 * Created by bittrich on 2/13/17.
 */
public class S03_RunKinkFinder {
    public static void main(String[] args) {
        runKinkFinder();
//        cleanDirectory();
    }

    /**
     * Screen for empty  result files and delete them. fix: actually do not delete them, so we can distinguish between
     * missing and failed computations.
     */
    private static void cleanDirectory() {
        Constants.list(Paths.get(Constants.KINK_FINDER_RESULT_PATH))
                .filter(path -> Constants.lines(path).count() == 1)
                .peek(path -> System.out.println("dropping empty file: " + path.toFile().getName()))
                .forEach(Constants::delete);
    }

    private static void runKinkFinder() {
        Path kinkFinderTmpPath = Paths.get(Constants.KINK_FINDER_TMP_PATH);
        Path kinkFinderResultPath = Paths.get(Constants.KINK_FINDER_RESULT_PATH);
        // based on all membrane protein structures
//        Stream.of("4NEF")
        Constants.lines(Paths.get(Constants.PDBTM_ALL_LIST))
                .map(line -> line.substring(0, 4))
                .distinct()
                .filter(pdbId -> !kinkFinderResultPath.resolve(pdbId + ".kinks").toFile().exists())
                .forEach(pdbId -> {
                    try {
                        System.out.println(pdbId);
                        // download file
                        URL url = new URL("https://files.rcsb.org/view/" + pdbId + ".pdb");
                        ReadableByteChannel rbc = Channels.newChannel(url.openStream());
                        String pdbFileName = kinkFinderTmpPath.resolve(pdbId + ".pdb").toFile().getAbsolutePath();
                        FileOutputStream fos = new FileOutputStream(pdbFileName);
                        fos.getChannel().transferFrom(rbc, 0, Long.MAX_VALUE);
                        fos.close();
                        rbc.close();

                        // run kink finder
                        ProcessBuilder processBuilder = new ProcessBuilder("python2",
                                Constants.KINK_FINDER_SCRIPT_PATH,
                                "-f",
                                pdbFileName,
                                "-o",
                                Constants.KINK_FINDER_TMP_PATH,
                                "-d");
                        processBuilder.inheritIO();
                        processBuilder.start().waitFor();
                        Files.copy(kinkFinderTmpPath.resolve("kinks.csv"), kinkFinderResultPath.resolve(pdbId + ".kinks"));
                    } catch (IOException e) {
                        throw new UncheckedIOException(e);
                    } catch (InterruptedException e) {
                        e.printStackTrace();
                    }
                });

        // based on non-redundant data set
//        Constants.list(Paths.get(Constants.STRUCTURE_PATH))
//                .filter(path -> !kinkFinderResultPath.resolve(path.toFile().getName().split("\\.")[0] + ".kinks").toFile().exists())
//                .map(Path::toFile)
//                .forEach(pdbPath -> {
//                    try {
//                        String pdbId = pdbPath.getName().split("\\.")[0];
//                        System.out.println(pdbId);
//                        ProcessBuilder processBuilder = new ProcessBuilder("python2",
//                                Constants.KINK_FINDER_SCRIPT_PATH,
//                                "-f",
//                                pdbPath.getAbsolutePath(),
//                                "-o",
//                                Constants.KINK_FINDER_TMP_PATH,
//                                "-d");
//                        processBuilder.inheritIO();
//                        processBuilder.start().waitFor();
//                        Files.copy(kinkFinderTmpPath.resolve("kinks.csv"), kinkFinderResultPath.resolve(pdbId + ".kinks"));
//                    } catch (IOException e) {
//                        throw new UncheckedIOException(e);
//                    } catch (InterruptedException e) {
//                        e.printStackTrace();
//                    }
//                });
    }
}
