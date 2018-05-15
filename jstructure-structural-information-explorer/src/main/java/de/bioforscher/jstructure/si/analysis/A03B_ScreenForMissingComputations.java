package de.bioforscher.jstructure.si.analysis;

import de.bioforscher.jstructure.mathematics.SetOperations;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class A03B_ScreenForMissingComputations {
    public static void main(String[] args) throws IOException {
        List<String> ids = Files.lines(Paths.get("/home/bittrich/git/phd_sb_repo/data/start2fold/pancsa-si.list"))
                .map(line -> line.split(";")[0])
                .collect(Collectors.toList());

        List<String> strategies = Stream.of(A03_ReconstructByVariousStrategy.ReconstructionStrategyDefinition.values())
                .map(A03_ReconstructByVariousStrategy.ReconstructionStrategyDefinition::getReconstructionStrategy)
                .map(A03_ReconstructByVariousStrategy.ReconstructionStrategy::getClass)
                .map(Class::getSimpleName)
                .collect(Collectors.toList());

        List<String> all = SetOperations.cartesianProductOf(ids, strategies)
                .map(pair -> pair.getLeft() + "," + pair.getRight())
                .collect(Collectors.toList());

        List<String> processed = Files.lines(Paths.get("/home/bittrich/git/phd_sb_repo/data/si/statistics/reconstruction-strategy.csv"))
                .filter(line -> line.startsWith("STF"))
                .map(line -> line.split(","))
                .map(split -> split[0] + "," + split[1])
                .distinct()
                .collect(Collectors.toList());

        List<String> missing = all.stream()
                .filter(combination -> !processed.contains(combination))
                .collect(Collectors.toList());

        missing.forEach(System.out::println);
    }
}
