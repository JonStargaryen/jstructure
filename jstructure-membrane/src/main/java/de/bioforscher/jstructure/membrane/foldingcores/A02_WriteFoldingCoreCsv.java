package de.bioforscher.jstructure.membrane.foldingcores;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.feature.energyprofile.EgorAgreement;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfile;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfileCalculator;
import de.bioforscher.jstructure.feature.interactions.*;
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
import java.util.Collection;
import java.util.Optional;
import java.util.stream.Collectors;

/**
 * For each residue in the dataset label:
 * <ul>
 *     <li>pdb-id</li>
 *     <li>chain-id</li>
 *     <li>res-id</li>
 *     <li>aa</li>
 *     <li>sse3</li>
 *     <li>sse9</li>
 *     <li>sseSize</li>
 *     <li># halogen</li>
 *     <li># hydrogen</li>
 *     <li># hydrophobic</li>
 *     <li># metal</li>
 *     <li># picat</li>
 *     <li># pistack</li>
 *     <li># salt</li>
 *     <li># water</li> //TODO could use shells, i.e. how many interactions are in the direct surroundings of the amino acid
 *     <li># total</li>
 *     <li># non-local hydrogen</li>
 *     <li># non-local hydrophobic</li>
 *     <li># non-local backbone</li>
 *     <li># non-local total</li>
 *     everything also normalized by max in the structure
 *     <!--li>sum PLIP energy</li--> //TODO could use rework or actual interaction energies
 *     <li>energy profile value</li>
 *     li>egor prediction</li>
 *     <li>DynaMine score</li>
 *     <li>early/late folding</li>
 * </ul>
 */
public class A02_WriteFoldingCoreCsv {
    private static final Logger logger = LoggerFactory.getLogger(A02_WriteFoldingCoreCsv.class);
    private static final DictionaryOfProteinSecondaryStructure DICTIONARY_OF_PROTEIN_SECONDARY_STRUCTURE = new DictionaryOfProteinSecondaryStructure();
    private static final EnergyProfileCalculator ENERGY_PROFILE_CALCULATOR = new EnergyProfileCalculator();
    private static final PLIPIntraMolecularAnnotator PLIP_INTRA_MOLECULAR_ANNOTATOR = new PLIPIntraMolecularAnnotator();
    private static final DynaMineBridge DYNA_MINE_BRIDGE = new DynaMineBridge();

    public static void main(String[] args) {
        String output = MembraneConstants.list(MembraneConstants.FOLDING_CORES_DIRECTORY.resolve("start2fold"))
                .map(A02_WriteFoldingCoreCsv::handleFile)
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.joining(System.lineSeparator(),
                        "pdb,chain,res,aa,sse3,sse9,sseSize," +
                                "freq_halogen,freq_hydrogen,freq_hydrophobic,freq_metal,freq_picat,freq_pistack,freq_salt,freq_water,freq_bb,freq_total," +
                                "freq_nl_hydrogen,freq_nl_hydrophobic,freq_nl_bb,freq_nl_total," +
                                "halogen,hydrogen,hydrophobic,metal,picat,pistack,salt,water,bb,total," +
                                "nl_hydrogen,nl_hydrophobic,nl_bb,nl_total," +
                                "energy,egor," +
                                "dynamine," +
                                "folds" + System.lineSeparator(),
                        ""));

        MembraneConstants.write(MembraneConstants.FOLDING_CORES_DIRECTORY.resolve("foldingcores.csv"),
                output);
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

                        PLIPInteractionContainer nonLocalInteractions = new PLIPInteractionContainer(null,
                                plipInteractionContainer
                                        .getInteractions()
                                        .stream()
                                        // interactions have to be non-local
                                        .filter(inter -> Math.abs(inter.getPartner1().getResidueIdentifier().getResidueNumber() - inter.getPartner2().getResidueIdentifier().getResidueNumber()) > 6)
                                        .collect(Collectors.toList()));
                        GenericSecondaryStructure sse = aminoAcid.getFeature(GenericSecondaryStructure.class);

                        double hal = determineMaxCount(chain, plipInteractionContainer, HalogenBond.class);
                        double hydrogen = determineMaxCount(chain, plipInteractionContainer, HydrogenBond.class);
                        double phobic = determineMaxCount(chain, plipInteractionContainer, HydrophobicInteraction.class);
                        double metal = determineMaxCount(chain, plipInteractionContainer, MetalComplex.class);
                        double picat = determineMaxCount(chain, plipInteractionContainer, PiCationInteraction.class);
                        double pistack = determineMaxCount(chain, plipInteractionContainer, PiStacking.class);
                        double salt = determineMaxCount(chain, plipInteractionContainer, SaltBridge.class);
                        double water = determineMaxCount(chain, plipInteractionContainer, WaterBridge.class);
                        double bb = determineMaxBackboneCount(chain, plipInteractionContainer);
                        double total = determineMaxCount(chain, plipInteractionContainer, PLIPInteraction.class);
                        double nl_hydro = determineMaxCount(chain, nonLocalInteractions, HydrogenBond.class);
                        double nl_phob = determineMaxCount(chain, nonLocalInteractions, HydrophobicInteraction.class);
                        double nl_bb = determineMaxBackboneCount(chain, nonLocalInteractions);
                        double nl_total = determineMaxCount(chain, nonLocalInteractions, PLIPInteraction.class);

                        return pdbId + "," +
                            chainId + "," +
                            aminoAcid.getResidueIdentifier() + "," +
                            aminoAcid.getOneLetterCode() + "," +
                            sse.getSecondaryStructure().getReducedRepresentation() + "," +
                            sse.getSecondaryStructure().getOneLetterRepresentation() + "," +
                            sse.getSurroundingSecondaryStructureElement(aminoAcid).getSize() + "," +
                            StandardFormat.format(plipInteractionContainer.getHalogenBonds().size() / hal) + "," +
                            StandardFormat.format(plipInteractionContainer.getHydrogenBonds().size() / hydrogen) + "," +
                            StandardFormat.format(plipInteractionContainer.getHydrophobicInteractions().size() / phobic) + "," +
                            StandardFormat.format(plipInteractionContainer.getMetalComplexes().size() / metal) + "," +
                            StandardFormat.format(plipInteractionContainer.getPiCationInteractions().size() / picat) + "," +
                            StandardFormat.format(plipInteractionContainer.getPiStackings().size() / pistack) + "," +
                            StandardFormat.format(plipInteractionContainer.getSaltBridges().size() / salt) + "," +
                            StandardFormat.format(plipInteractionContainer.getWaterBridges().size() / water) + "," +
                            StandardFormat.format(plipInteractionContainer.getBackboneInteractions().size() / bb) + "," +
                            StandardFormat.format(plipInteractionContainer.getInteractions().size() / total) + "," +
                            StandardFormat.format(nonLocalInteractions.getHydrogenBonds().size() / nl_hydro) + "," +
                            StandardFormat.format(nonLocalInteractions.getHydrophobicInteractions().size() / nl_phob) + "," +
                            StandardFormat.format(nonLocalInteractions.getBackboneInteractions().size() / nl_bb) + "," +
                            StandardFormat.format(nonLocalInteractions.getInteractions().size() / nl_total)+ "," +
                            plipInteractionContainer.getHalogenBonds().size() + "," +
                            plipInteractionContainer.getHydrogenBonds().size() + "," +
                            plipInteractionContainer.getHydrophobicInteractions().size() + "," +
                            plipInteractionContainer.getMetalComplexes().size() + "," +
                            plipInteractionContainer.getPiCationInteractions().size() + "," +
                            plipInteractionContainer.getPiStackings().size() + "," +
                            plipInteractionContainer.getSaltBridges().size() + "," +
                            plipInteractionContainer.getWaterBridges().size() + "," +
                            plipInteractionContainer.getBackboneInteractions().size() + "," +
                            plipInteractionContainer.getInteractions().size() + "," +
                            nonLocalInteractions.getHydrogenBonds().size() + "," +
                            nonLocalInteractions.getHydrophobicInteractions().size() + "," +
                            nonLocalInteractions.getBackboneInteractions().size() + "," +
                            nonLocalInteractions.getInteractions().size() + "," +
                            StandardFormat.format(aminoAcid.getFeature(EnergyProfile.class).getSolvationEnergy()) + "," +
                            StandardFormat.format(aminoAcid.getFeature(EgorAgreement.class).getEgorPrediction()) + "," +
                            StandardFormat.format(aminoAcid.getFeature(BackboneRigidity.class).getBackboneRigidity()) + "," + 
                            earlyLate;
                    })
                    .collect(Collectors.joining(System.lineSeparator())));
        } catch (Exception e) {
            logger.warn("computation for {} failed",
                    pdbId,
                    e);
            return Optional.empty();
        }
    }

    private static double determineMaxCount(Chain chain, PLIPInteractionContainer container, Class<? extends PLIPInteraction> interaction) {
        double count = chain.aminoAcids()
                .map(container::getInteractionsFor)
                .map(PLIPInteractionContainer::getInteractions)
                .flatMap(Collection::stream)
                .filter(interaction::isInstance)
                .count();
        return count > 0 ? count : 1;
    }

    private static double determineMaxBackboneCount(Chain chain, PLIPInteractionContainer container) {
        double count = chain.aminoAcids()
                .map(container::getInteractionsFor)
                .map(PLIPInteractionContainer::getBackboneInteractions)
                .mapToLong(Collection::size)
                .sum();
        return count > 0 ? count : 1;
    }
}
