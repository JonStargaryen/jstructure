package de.bioforscher.jstructure.membrane.foldingcores;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceArea;
import de.bioforscher.jstructure.feature.energyprofile.EgorAgreement;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfile;
import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.membrane.MembraneConstants;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Optional;
import java.util.stream.Collectors;

public class A03_WriteFoldingCoreCsvWithConventionalContactDefinition {
    private static final Logger logger = LoggerFactory.getLogger(A03_WriteFoldingCoreCsvWithConventionalContactDefinition.class);

    public static void main(String[] args) {
        String output = MembraneConstants.list(MembraneConstants.FOLDING_CORES_DIRECTORY.resolve("start2fold"))
                .map(A03_WriteFoldingCoreCsvWithConventionalContactDefinition::handleFile)
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.joining(System.lineSeparator(),
                        "pdb,chain,res,aa,sse3,sse9,sseSize," +
                                "contacts,nl_contacts," +
                                "energy,egor," +
                                "asa," +
                                "eccount,cumstrength,ecstrength,conservation," +
                                "folds" + System.lineSeparator(),
                        ""));

        MembraneConstants.write(MembraneConstants.FOLDING_CORES_DIRECTORY.resolve("foldingcores-conv.csv"),
                output);
    }

    private static final double BETA_THRESHOLD = 8.0;

    private static Optional<String> handleFile(Path path) {
        try {
            String id = Files.lines(path)
                    .filter(line -> line.startsWith("#pdb:"))
                    .findFirst()
                    .get()
                    .split(": ")[1];
            String pdbId = id.split("_")[0];
            String chainId = id.split("_")[1];
            logger.info("processing {} chain {}",
                    pdbId,
                    chainId);

            Structure structure = StructureParser.source(pdbId)
                    .parse();
            Chain chain = structure.select()
                    .chainId(chainId)
                    .asChain();

            // parse coupling scores

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
                        EvolutionaryCouplingParser.parseHotSpotFile(chain,
                                MembraneConstants.FOLDING_CORES_DIRECTORY.resolve("couplings").resolve(pdbId + "_A_hs.html"));

                        EvolutionaryCouplingParser.HotSpotScoring hotSpotScoring = aminoAcid.getFeature(EvolutionaryCouplingParser.HotSpotScoring.class);

                        return pdbId + "," +
                                chainId + "," +
                                aminoAcid.getResidueIdentifier() + "," +
                                aminoAcid.getOneLetterCode() + "," +
                                sse.getSecondaryStructure().getReducedRepresentation() + "," +
                                sse.getSecondaryStructure().getOneLetterRepresentation() + "," +
                                sse.getSurroundingSecondaryStructureElement(aminoAcid).getSize() + "," +
                                interactions + "," +
                                nonLocalInteractions + "," +
                                StandardFormat.format(aminoAcid.getFeature(EnergyProfile.class).getSolvationEnergy()) + "," +
                                StandardFormat.format(aminoAcid.getFeature(EgorAgreement.class).getEgorPrediction()) + "," +
                                StandardFormat.format(aminoAcid.getFeature(AccessibleSurfaceArea.class).getRelativeAccessibleSurfaceArea()) + "," +
                                hotSpotScoring.getEcCount() + "," +
                                StandardFormat.format(hotSpotScoring.getCumStrength()) + "," +
                                StandardFormat.format(hotSpotScoring.getEcStrength()) + "," +
                                hotSpotScoring.getConservation() + "," +
                                earlyLate;
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
}
