package de.bioforscher.jstructure.membrane.foldingcores;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceArea;
import de.bioforscher.jstructure.feature.energyprofile.EgorAgreement;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfile;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.mathematics.LinearAlgebra;
import de.bioforscher.jstructure.membrane.MembraneConstants;
import de.bioforscher.jstructure.model.structure.*;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class A03_WriteFoldingCoreCsvWithConventionalContactDefinition {
    private static final Logger logger = LoggerFactory.getLogger(A03_WriteFoldingCoreCsvWithConventionalContactDefinition.class);

    public static void main(String[] args) {
        String output = MembraneConstants.list(MembraneConstants.FOLDING_CORES_DIRECTORY.resolve("start2fold"))
                .map(A03_WriteFoldingCoreCsvWithConventionalContactDefinition::handleFile)
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.joining(System.lineSeparator(),
                        "pdb,chain,res,aa,sse3,sse9,sseSize," +
                                "rasa_core,ep_core," +
                                // conventional contact definition
                                "contacts,nl_contacts," +
                                // PLIP contact definition
                                "plip_halogen,plip_hydrogen,plip_hydrophobic,plip_metal,plip_picat,plip_pistack,plip_salt,plip_water,plip_bb,plip_total," +
                                "plip_nl_hydrogen,plip_nl_hydrophobic,plip_nl_bb,plip_nl_total," +
                                "energy,egor," +
                                "asa," +
                                "eccount,cumstrength,ecstrength,conservation," +
                                "folds" + System.lineSeparator(),
                        ""));

        MembraneConstants.write(MembraneConstants.FOLDING_CORES_DIRECTORY.resolve("foldingcores.csv"),
                output);
    }

    private static final double BETA_THRESHOLD = 8.0;

    private static Optional<String> handleFile(Path path) {
        try {
            String id;
            try(Stream<String> lines = Files.lines(path)) {
                id = lines.filter(line -> line.startsWith("#pdb:"))
                        .findFirst()
                        .get()
                        .split(": ")[1];
            }
            String pdbId = id.split("_")[0];
            String chainId = id.split("_")[1];
            logger.info("processing {} chain {}",
                    pdbId,
                    chainId);

            Structure structure = StructureParser.source(pdbId).parse();
            Chain chain = structure.select()
                    .chainId(chainId)
                    .asChain();

            EvolutionaryCouplingParser.parseHotSpotFile(chain,
                    MembraneConstants.FOLDING_CORES_DIRECTORY.resolve("couplings").resolve(pdbId + "_A_hs.html"));

            return Optional.of(chain.aminoAcids()
                    .map(aminoAcid -> {
                        try {
                            String earlyLate;
                            try (Stream<String> lines = Files.lines(path)) {
                                earlyLate = lines.filter(line -> !line.startsWith("#"))
                                        .filter(line -> line.endsWith("EARLY"))
                                        .map(line -> line.split(";")[0])
                                        .filter(resNum -> aminoAcid.getResidueIdentifier().toString().equals(resNum))
                                        .findFirst()
                                        .map(resNum -> "early")
                                        .orElse("late");
                            }

                            int interactions = (int) chain.aminoAcids()
                                    // ignore amino acid itself
                                    .filter(aa -> !aminoAcid.equals(aa))
                                    .filter(aa -> getBetaCarbon(aminoAcid).calculate().distance(getBetaCarbon(aa)) < BETA_THRESHOLD)
                                    .count();

                            int nonLocalInteractions = (int) chain.aminoAcids()
                                    // interactions have to be non-local
                                    .filter(aa -> Math.abs(aminoAcid.getResidueIdentifier().getResidueNumber() - aa.getResidueIdentifier().getResidueNumber()) > 6)
                                    .filter(aa -> getBetaCarbon(aminoAcid).calculate().distance(getBetaCarbon(aa)) < BETA_THRESHOLD)
                                    .count();

                            GenericSecondaryStructure sse = aminoAcid.getFeature(GenericSecondaryStructure.class);

                            EvolutionaryCouplingParser.HotSpotScoring hotSpotScoring = aminoAcid.getFeature(EvolutionaryCouplingParser.HotSpotScoring.class);

                            PLIPInteractionContainer plipInteractionContainer = aminoAcid.getFeature(PLIPInteractionContainer.class);
                            PLIPInteractionContainer nonLocalPlipInteractionContainer = new PLIPInteractionContainer(null,
                                    plipInteractionContainer
                                            .getInteractions()
                                            .stream()
                                            // interactions have to be non-local
                                            .filter(inter -> Math.abs(inter.getPartner1().getResidueIdentifier().getResidueNumber() - inter.getPartner2().getResidueIdentifier().getResidueNumber()) > 6)
                                            .collect(Collectors.toList()));

                            return pdbId + "," +
                                    chainId + "," +
                                    aminoAcid.getResidueIdentifier() + "," +
                                    aminoAcid.getOneLetterCode() + "," +
                                    sse.getSecondaryStructure().getReducedRepresentation() + "," +
                                    sse.getSecondaryStructure().getOneLetterRepresentation() + "," +
                                    sse.getSurroundingSecondaryStructureElement(aminoAcid).getSize() + "," +
                                    applyRasaInsideOutsideCriterion(aminoAcid) + "," +
                                    applyEnergyProfileInsideOutsideCriterion(aminoAcid, chain) + "," +
                                    interactions + "," +
                                    nonLocalInteractions + "," +
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
                                    nonLocalPlipInteractionContainer.getHydrogenBonds().size() + "," +
                                    nonLocalPlipInteractionContainer.getHydrophobicInteractions().size() + "," +
                                    nonLocalPlipInteractionContainer.getBackboneInteractions().size() + "," +
                                    nonLocalPlipInteractionContainer.getInteractions().size() + "," +
                                    StandardFormat.format(aminoAcid.getFeature(EnergyProfile.class).getSolvationEnergy()) + "," +
                                    StandardFormat.format(aminoAcid.getFeature(EgorAgreement.class).getEgorPrediction()) + "," +
                                    StandardFormat.format(aminoAcid.getFeature(AccessibleSurfaceArea.class).getRelativeAccessibleSurfaceArea()) + "," +
                                    hotSpotScoring.getEcCount() + "," +
                                    StandardFormat.format(hotSpotScoring.getCumStrength()) + "," +
                                    StandardFormat.format(hotSpotScoring.getEcStrength()) + "," +
                                    hotSpotScoring.getConservation() + "," +
                                    earlyLate;
                        } catch (IOException e) {
                            throw new UncheckedIOException(e);
                        }
                    })
                    .collect(Collectors.joining(System.lineSeparator())));
        } catch (Exception e) {
            logger.warn("computation for {} failed",
                    path,
                    e);
            return Optional.empty();
        }
    }

    private static Atom getBetaCarbon(AminoAcid aminoAcid) {
        return aminoAcid.atoms()
                .filter(atom -> atom.getName().equals("CB"))
                .findFirst()
                .orElse(aminoAcid.getCa());
    }

    /**
     * see Rost B, Sander C. Conservation and prediction of solvent accessibility in protein families. Proteins 1994 Nov;20(3):216-26.
     * @param aminoAcid
     * @return
     */
    private static boolean applyRasaInsideOutsideCriterion(AminoAcid aminoAcid) {
        return aminoAcid.getFeature(AccessibleSurfaceArea.class).getRelativeAccessibleSurfaceArea() < 0.16;
    }

    /**
     * see Frank Dressel's PhD Thesis
     * @param aminoAcid the considered amino acid
     * @param chain all amino acids of this chain
     * @return true if residue is inside according to definition, false if outside
     */
    private static boolean applyEnergyProfileInsideOutsideCriterion(AminoAcid aminoAcid, Chain chain) {
        Atom alphaCarbon = aminoAcid.getCa();
        Atom betaCarbon = getBetaCarbon(aminoAcid);
        double[] c = chain.aminoAcids()
                .filter(aa -> aa.getCa().calculate().distance(alphaCarbon) < 10)
                .map(AminoAcid::getCa)
                .collect(StructureCollectors.toCentroid());
               // | C_alpha - c | < 5
        return (Math.abs(LinearAlgebra.on(alphaCarbon).subtract(c).norm()) < 5) &&
                // ( C_alpha - C_beta ) ( C_alpha - c ) < 0
                ((LinearAlgebra.on(alphaCarbon).subtract(betaCarbon))
                        .dotProduct(LinearAlgebra.on(alphaCarbon).subtract(c)) < 0);
    }
}
