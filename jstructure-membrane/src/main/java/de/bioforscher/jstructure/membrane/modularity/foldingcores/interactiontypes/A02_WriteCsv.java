package de.bioforscher.jstructure.membrane.modularity.foldingcores.interactiontypes;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfile;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfileCalculator;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.interactions.PLIPIntraMolecularAnnotator;
import de.bioforscher.jstructure.feature.rigidity.BackboneRigidity;
import de.bioforscher.jstructure.feature.rigidity.DynaMineBridge;
import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.feature.sse.dssp.DictionaryOfProteinSecondaryStructure;
import de.bioforscher.jstructure.membrane.MembraneConstants;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Optional;
import java.util.stream.Collectors;

/**
 * For each residue in the dataset label:
 * <ul>
 *     <li>pdb-id</li>
 *     <li>chain-id</li>
 *     <li>res-id</li>
 *     <li>aa</li>
 *     <li>sse</li>
 *     <li>early/late folding</li>
 *     <li># halogen</li>
 *     <li># hydrogen</li>
 *     <li># hydrophobic</li>
 *     <li># metal</li>
 *     <li># picat</li>
 *     <li># pistack</li>
 *     <li># salt</li>
 *     <li># water</li>
 *     <li># total</li>
 *     <li>sum PLIP energy</li>
 *     <li>energy profile value</li>
 *     <li>DynaMine score</li>
 * </ul>
 */
public class A02_WriteCsv {
    private static final Logger logger = LoggerFactory.getLogger(A02_WriteCsv.class);
    private static final DictionaryOfProteinSecondaryStructure DICTIONARY_OF_PROTEIN_SECONDARY_STRUCTURE = new DictionaryOfProteinSecondaryStructure();
    private static final EnergyProfileCalculator ENERGY_PROFILE_CALCULATOR = new EnergyProfileCalculator();
    private static final PLIPIntraMolecularAnnotator PLIP_INTRA_MOLECULAR_ANNOTATOR = new PLIPIntraMolecularAnnotator();
    private static final DynaMineBridge DYNA_MINE_BRIDGE = new DynaMineBridge();

    public static void main(String[] args) {
        String output = MembraneConstants.list(MembraneConstants.FOLDING_CORES_DIRECTORY.resolve("start2fold"))
                .map(A02_WriteCsv::handleFile)
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.joining(System.lineSeparator(),
                        "pdb;chain;res;aa;sse;folds;halogen;hydrogen;hydrophobic;metal;picat;pistack;salt;water;total;energy;dynamine" + System.lineSeparator(),
                        ""));

        System.out.println(output);
    }

    private static Optional<String> handleFile(Path path) {
        String id = MembraneConstants.lines(path)
                .filter(line -> line.startsWith("#pdb:"))
                .findFirst()
                .get()
                .split(": ")[1];
        String pdbId = id.split("_")[0];
        String chainId = id.split("_")[1];
        logger.info("processing {} chain {}",
                pdbId,
                chainId);

        try {
            Structure structure = StructureParser.source(pdbId)
                    .minimalParsing(true)
                    .parse();
            Chain chain = structure.select()
                    .chainId(chainId)
                    .asChain();

            logger.info("annotating secondary structure elements");
            DICTIONARY_OF_PROTEIN_SECONDARY_STRUCTURE.process(structure);
            logger.info("computing energy profile");
            ENERGY_PROFILE_CALCULATOR.process(structure);
            logger.info("annotating intra-molecular interactions");
            PLIP_INTRA_MOLECULAR_ANNOTATOR.process(structure);
            logger.info("parsing dynamine data");
            DYNA_MINE_BRIDGE.process(chain, Files.lines(path.getParent()
                    .getParent()
                    .resolve("dynamine")
                    .resolve(pdbId + "_" + chainId + "_" + path.toFile().getName().split("\\.")[0] + "_" + path.toFile().getName().split("\\.")[1] + "_backbone.pred"))
                    .collect(Collectors.joining(System.lineSeparator())));

            return Optional.of(chain.aminoAcids()
                    .map(aminoAcid -> {
                        String earlyLate = MembraneConstants.lines(path)
                                .filter(line -> !line.startsWith("#"))
                                .filter(line -> line.endsWith("EARLY"))
                                .map(line -> line.split(";")[0])
                                .filter(resNum -> aminoAcid.getResidueIdentifier().toString().equals(resNum))
                                .findFirst()
                                .map(resNum -> "early")
                                .orElse("late");
                        PLIPInteractionContainer plipInteractionContainer = aminoAcid.getFeature(PLIPInteractionContainer.class);
                        return pdbId + ";" +
                            chainId + ";" +
                            aminoAcid.getResidueIdentifier() + ";" +
                            aminoAcid.getOneLetterCode() + ";" +
                            aminoAcid.getFeature(GenericSecondaryStructure.class).getSecondaryStructure().getOneLetterRepresentation() + ";" +
                            earlyLate + ";" +
                            plipInteractionContainer.getHalogenBonds().size() + ";" +
                            plipInteractionContainer.getHydrogenBonds().size() + ";" +
                            plipInteractionContainer.getHydrophobicInteractions().size() + ";" +
                            plipInteractionContainer.getMetalComplexes().size() + ";" +
                            plipInteractionContainer.getPiCationInteractions().size() + ";" +
                            plipInteractionContainer.getPiStackings().size() + ";" +
                            plipInteractionContainer.getSaltBridges().size() + ";" +
                            plipInteractionContainer.getWaterBridges().size() + ";" +
                            plipInteractionContainer.getInteractions().size() + ";" +
                            StandardFormat.format(aminoAcid.getFeature(EnergyProfile.class).getSolvationEnergy()) + ";" +
                            aminoAcid.getFeature(BackboneRigidity.class).getBackboneRigidity();
                    })
                    .collect(Collectors.joining(System.lineSeparator())));
        } catch (Exception e) {
            logger.warn("computation for {} failed",
                    pdbId,
                    e);
            return Optional.empty();
        }
    }
}
