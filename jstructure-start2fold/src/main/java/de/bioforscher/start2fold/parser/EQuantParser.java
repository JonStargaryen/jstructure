package de.bioforscher.start2fold.parser;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.start2fold.model.EQuantScore;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.stream.Stream;

public class EQuantParser {
    public static void parseEQuantFile(Chain chain,
                                       Path equantPath) {
        try {
            String chainId = chain.getChainIdentifier().getChainId();
            try (Stream<String> lines = Files.lines(equantPath)) {
                // skip header
                lines.filter(line -> !line.startsWith("chain"))
                        .filter(line -> line.startsWith(chainId))
                        .forEach(line -> {
                            String[] split = line.split("\\s+");
                            int residueNumber = Integer.valueOf(split[1]);
                            double evaluation = Double.valueOf(split[4]);
                            AminoAcid aminoAcid = chain.select()
                                    .residueNumber(residueNumber)
                                    .asAminoAcid();
                            aminoAcid.getFeatureContainer().addFeature(new EQuantScore(evaluation));
                        });
            }
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }
}
