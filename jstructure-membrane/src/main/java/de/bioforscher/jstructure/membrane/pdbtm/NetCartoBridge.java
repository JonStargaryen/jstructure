package de.bioforscher.jstructure.membrane.pdbtm;

import de.bioforscher.jstructure.membrane.MembraneConstants;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.annotation.concurrent.NotThreadSafe;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * netcarto_cl testnetwork.dat 1111 5.0 0.4 0.96 5
 *             input           seed T_i i   c    N_r
 *             the input file
 *                             for stochastic process
 *                                  initial temperature
 *                                      iteration factor: high values more precise, but slower
 *                                          cooling factor: again high values more precise, but slower
 *                                               number of randomizations: evaluates significance, can be 0
 *
 * Even though this class is stateless, NetCarlo creates temporary files, thus, creating a state on filesystem level.
 */
@NotThreadSafe
public class NetCartoBridge {
    private static final Logger logger = LoggerFactory.getLogger(NetCartoBridge.class);
    private static final List<String> FILESNAMES = Stream.of("modules.clu",
            "modules.dat",
            "network.net",
            "node_prop.dat",
            "randomized_mod.dat",
            "roles.clu",
            "roles.dat").collect(Collectors.toList());
    private final Path workingDirectory;

    NetCartoBridge(Path workingDirectory) {
        this.workingDirectory = workingDirectory;

        MembraneConstants.list(workingDirectory)
                // skip executable etc
                .filter(path -> path.toFile().getName().endsWith(".dat"))
                //TODO skip naive networks for now
                .filter(path -> path.toFile().getName().endsWith("_plip.dat"))
                .forEach(this::processNetworkFile);
    }

    private void processNetworkFile(Path path) {
        String name = path.toFile().getName();
        logger.info("processing {}", name);
        String id = name.split("\\.")[0];

        try {
            String[] command = new String[] { workingDirectory.resolve("netcarto_cl.exe").toFile().getAbsolutePath(),
                    workingDirectory.resolve(name).toFile().getAbsolutePath(),
                    "1111",
                    "5.0",
                    "0.4",
                    "0.96",
                    "5"
            };
            logger.info("{}", Stream.of(command).collect(Collectors.joining(" ")));
            new ProcessBuilder(command)
                    .inheritIO()
                    .start()
                    .waitFor();

            for(String outputFile : FILESNAMES) {
                try {
                    Path outputPath = Paths.get("C:/Users/S/git/jstructure/").resolve(outputFile);
                    Path renamedOutputPath = workingDirectory.resolve(id + "." + outputFile);
                    logger.info("moving result file {} => {}",
                            outputPath.toFile().getName(),
                            renamedOutputPath.toFile().getName());
                    Files.move(outputPath, renamedOutputPath);
                } catch (IOException e){
                    logger.warn("failed to rename temp file {}",
                            outputFile,
                            e);
                }
            }
        } catch (Exception e) {
            logger.warn("failed to compute modularity for {}",
                    name,
                    e);
        }
    }
}
