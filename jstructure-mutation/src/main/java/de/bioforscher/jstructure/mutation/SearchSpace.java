package de.bioforscher.jstructure.mutation;

import de.bioforscher.jstructure.feature.mapping.SiftsMappingAnnotator;
import de.bioforscher.jstructure.model.structure.identifier.ChainIdentifier;
import de.bioforscher.jstructure.model.structure.identifier.ProteinIdentifier;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Defines a search space to explore.
 * Created by bittrich on 7/6/17.
 */
public class SearchSpace {
    public static List<ProteinIdentifier> fullPdb() {
        try {
            return Files.walk(Paths.get("/home/bittrich/pdb/"))
                    .filter(path -> !Files.isDirectory(path))
                    .map(Path::toFile)
                    .map(File::getName)
                    .map(name -> name.substring(3, 7))
                    .map(ProteinIdentifier::createFromPdbId)
                    .collect(Collectors.toList());
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    public static List<ChainIdentifier> pdbByPfam(String pfamId) {
        return SiftsMappingAnnotator.getLinesForPfamId(pfamId)
                .map(split -> ChainIdentifier.createFromChainId(ProteinIdentifier.createFromPdbId(split[0]), split[1]))
                .collect(Collectors.toList());
    }

    public static List<ChainIdentifier> pdbByUniProt(String uniProtId) {
        return SiftsMappingAnnotator.getLinesForUniProtId(uniProtId)
                .map(split -> ChainIdentifier.createFromChainId(ProteinIdentifier.createFromPdbId(split[0]), split[1]))
                .collect(Collectors.toList());
    }
}
