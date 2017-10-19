package de.bioforscher.jstructure.membrane.modularity.foldingcores;

import de.bioforscher.jstructure.membrane.MembraneConstants;

import java.nio.file.Path;

/**
 * Based on: http://www.sciencedirect.com/science/article/pii/S0006349515048122?via%3Dihub#mmc1
 */
public class A01_CheckSanity {
    public static void main(String[] args) {
        MembraneConstants.lines(MembraneConstants.FOLDING_CORES_DIRECTORY.resolve("sanity.txt"))
                .filter(line -> !line.startsWith("#"))
                .forEach(A01_CheckSanity::checkSanity);
    }

    private static void checkSanity(String line) {
//        System.out.println(line);
        String pdbId = line.split(",")[0];
        int numberOfEarlyFoldingResidues = Integer.valueOf(line.split(",")[1]);

        Path start2FoldFilePath = MembraneConstants.list(MembraneConstants.FOLDING_CORES_DIRECTORY.resolve("start2fold"))
                .filter(path -> MembraneConstants.lines(path)
                        .anyMatch(l -> l.contains(pdbId)))
                .findFirst()
                .get();

        long obs = MembraneConstants.lines(start2FoldFilePath)
                .filter(l -> !l.startsWith("#"))
                .filter(l -> l.endsWith("EARLY"))
                .distinct()
                .count();

        if(numberOfEarlyFoldingResidues != obs) {
            System.err.println(pdbId + " - " + start2FoldFilePath.toFile().getName() +  " - expected: " + numberOfEarlyFoldingResidues + " - found: " + obs);
        }
    }
}
