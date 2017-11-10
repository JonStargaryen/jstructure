package de.bioforscher.jstructure.membrane.weblogos;

import de.bioforscher.jstructure.membrane.MembraneConstants;

import java.nio.file.Path;
import java.util.stream.Collectors;

public class A02_DownsizeSequences {
    public static void main(String[] args) {
        MembraneConstants.list(MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY
                .resolve("results")
                .resolve("weblogos"))
                .filter(path -> !path.toFile().getName().contains("small"))
                .forEach(originalPath -> {
                    String name = originalPath.toFile().getName();
                    Path outputPath = originalPath.getParent().resolve(name + "-small");
                    System.out.println(name);

                    MembraneConstants.list(originalPath)
                            .forEach(file -> {
                                String output = MembraneConstants.lines(file)
                                        .map(line -> line.substring(6, 15))
                                        .collect(Collectors.joining(System.lineSeparator()));
                                MembraneConstants.write(outputPath.resolve(file.toFile().getName()), output);
                            });
                });
    }
}
