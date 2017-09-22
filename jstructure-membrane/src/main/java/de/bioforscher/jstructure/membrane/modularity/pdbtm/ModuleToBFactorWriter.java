package de.bioforscher.jstructure.membrane.modularity.pdbtm;

import de.bioforscher.jstructure.membrane.MembraneConstants;
import de.bioforscher.jstructure.model.identifier.IdentifierFactory;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.nio.file.Path;
import java.util.List;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

public class ModuleToBFactorWriter {
    private static final Logger logger = LoggerFactory.getLogger(ModuleToBFactorWriter.class);
    private final Path pdbPath;
    private final Path networkPath;
    private final Path bfactorPath;

    public ModuleToBFactorWriter(Path datasetPath) {
        this.pdbPath = datasetPath.resolve("pdb");
        this.networkPath = datasetPath.resolve("network");
        this.bfactorPath = networkPath.resolve("bfactor");

        MembraneConstants.list(networkPath)
                .filter(path -> path.toFile().getName().endsWith(".modules.dat"))
                .forEach(this::handleFile);
    }

    private void handleFile(Path path) {
        String name = path.toFile().getName();
        logger.info("handling {}", name);
        String pdbId = name.split("_")[0];
        String chainId = name.split("_")[1];

        Structure chain = StructureParser.source(pdbPath.resolve(pdbId + ".pdb"))
                .minimalParsing(true)
                .parse()
                .select()
                .chainId(chainId)
                .asIsolatedStructure();

        // set bfactor to default case - i.e. no module assigned
        chain.atoms().forEach(atom -> atom.setBfactor(0));

        // fix 09/20/17: ignore modules < 10 nodes, use fraction to represent clusters (enables globally comparable coloring)
        int numberOfSignificantModules = (int) MembraneConstants.lines(path)
                .filter(line -> !line.startsWith("#"))
                .filter(line -> line.split("---")[1].trim().split("\\s+").length > 9)
                .count();
        int[] processedClusters = { 0 };

        MembraneConstants.lines(path)
                .filter(line -> !line.startsWith("#"))
                .forEach(line -> {
                    int moduleNumber = Integer.valueOf(line.split("\\s+")[0]);
                    List<ResidueIdentifier> residueIdentifiers = Pattern.compile("\\s+")
                            .splitAsStream(line.split("---")[1].trim())
                            .map(IdentifierFactory::createResidueIdentifier)
                            .collect(Collectors.toList());
                    logger.info("module '{}' contains nodes: {}",
                            moduleNumber,
                            residueIdentifiers);

                    if(residueIdentifiers.size() < 10) {
                        logger.info("skipping sparsely populated module");
                        return;
                    }

                    processedClusters[0]++;

                    for(ResidueIdentifier residueIdentifier : residueIdentifiers) {
                        Group group = chain.select()
                                .residueIdentifier(residueIdentifier)
                                .asGroup();
                        group.atoms().forEach(atom -> atom.setBfactor((float) processedClusters[0] / (float) numberOfSignificantModules));
                    }
                });

        // write manipulated structure
        Path outputPath = bfactorPath.resolve(name + ".pdb");
        logger.info("writing output file {}", outputPath);
        MembraneConstants.write(outputPath, chain.getPdbRepresentation());
    }

//    private void handleId(String id) {
//        logger.info("handling {}", id);
//        String pdbId = id.split("_")[0];
//        String chainId = id.split("_")[1];
//
//        Structure chain = StructureParser.source(pdbPath.resolve(pdbId + ".pdb"))
//                .minimalParsing(true)
//                .parse()
//                .select()
//                .chainId(chainId)
//                .asIsolatedStructure();
//
//        // set bfactor to default case - i.e. no module assigned
//        chain.atoms().forEach(atom -> atom.setBfactor(0));
//
//        MembraneConstants.lines(networkPath.resolve(id + "_plip.modules.dat"))
//                .filter(line -> !line.startsWith("#"))
//                .forEach(line -> {
//                    int moduleNumber = Integer.valueOf(line.split("\\s+")[0]);
//                    List<ResidueIdentifier> residueIdentifiers = Pattern.compile("\\s+")
//                            .splitAsStream(line.split("---")[1].trim())
//                            .map(IdentifierFactory::createResidueIdentifier)
//                            .collect(Collectors.toList());
//                    logger.info("module '{}' contains nodes: {}",
//                            moduleNumber,
//                            residueIdentifiers);
//
//                    for(ResidueIdentifier residueIdentifier : residueIdentifiers) {
//                        Group group = chain.select()
//                                .residueIdentifier(residueIdentifier)
//                                .asGroup();
//                        group.atoms().forEach(atom -> atom.setBfactor(moduleNumber));
//                    }
//                });
//
//        // write manipulated structure
//        Path outputPath = bfactorPath.resolve(id + ".pdb");
//        logger.info("writing output file {}", outputPath);
//        MembraneConstants.write(outputPath, chain.getPdbRepresentation());
//    }
}
