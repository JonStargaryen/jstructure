package de.bioforscher.jstructure.gmlvq;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceArea;
import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceAreaCalculator;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfile;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfileCalculator;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.interactions.PLIPIntraMolecularAnnotator;
import de.bioforscher.jstructure.feature.loopfraction.LoopFraction;
import de.bioforscher.jstructure.feature.loopfraction.LoopFractionCalculator;
import de.bioforscher.jstructure.feature.sse.dssp.DictionaryOfProteinSecondaryStructure;
import de.bioforscher.jstructure.model.identifier.IdentifierFactory;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.testutil.FileUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * Create the class-less arff file for the HDS dataset of the 3DSIG poster.
 */
public class A01_AnnotateFeatures {
    private static final Logger logger = LoggerFactory.getLogger(A01_AnnotateFeatures.class);
    private static final double P_VALUE_CUTOFF = 0.001;
    private static int counter;
    private static int binSize;
    private static final AccessibleSurfaceAreaCalculator ACCESSIBLE_SURFACE_AREA_CALCULATOR =
            new AccessibleSurfaceAreaCalculator();
    private static final EnergyProfileCalculator ENERGY_PROFILE_CALCULATOR =
            new EnergyProfileCalculator();
    private static final DictionaryOfProteinSecondaryStructure DICTIONARY_OF_PROTEIN_SECONDARY_STRUCTURE =
            new DictionaryOfProteinSecondaryStructure();
    private static final LoopFractionCalculator LOOP_FRACTION_CALCULATOR =
            new LoopFractionCalculator();
    private static final PLIPIntraMolecularAnnotator PLIP_INTRA_MOLECULAR_ANNOTATOR =
            new PLIPIntraMolecularAnnotator();

    public static void main(String[] args) throws IOException {
        Path summaryPath = FileUtils.GIT_DIRECTORY.resolve("gmlvq_main")
                .resolve("data")
                .resolve("csa_new")
                .resolve("motifs_BLAST_10e80_matches")
                .resolve("1a0j_1")
                .resolve("summary.csv");

        Map<String, List<String>> proteinMap = Files.lines(summaryPath)
                .filter(line -> !line.startsWith("match"))
                .filter(line -> Double.valueOf(line.split(",")[2]) < P_VALUE_CUTOFF)
                .collect(Collectors.groupingBy(line -> line.split("_")[0]));

        counter = 1;
        binSize = proteinMap.size();

        String output = proteinMap.entrySet()
                .stream()
                .map(A01_AnnotateFeatures::handleBin)
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.joining(System.lineSeparator()));

        Files.write(summaryPath.getParent().resolve("summary.arff"), output.getBytes());
    }

    private static Optional<String> handleBin(Map.Entry<String, List<String>> entry) {
        String pdbId = entry.getKey();
        try {
            logger.info("annotating {} - {} / {}",
                    pdbId,
                    counter,
                    binSize);
            Structure structure = StructureParser.source(pdbId)
                    .minimalParsing(true)
                    .parse();

            PLIP_INTRA_MOLECULAR_ANNOTATOR.process(structure);
            DICTIONARY_OF_PROTEIN_SECONDARY_STRUCTURE.process(structure);
            LOOP_FRACTION_CALCULATOR.process(structure);
            ACCESSIBLE_SURFACE_AREA_CALCULATOR.process(structure);
            ENERGY_PROFILE_CALCULATOR.process(structure);

            return Optional.of(entry.getValue()
                    .stream()
                    .map(hit -> handleHit(structure, hit))
                    .filter(Optional::isPresent)
                    .map(Optional::get)
                    .collect(Collectors.joining(System.lineSeparator())));
        } catch (Exception e) {
            e.printStackTrace();
            return Optional.empty();
        } finally {
            counter++;
        }
    }

    private static Optional<String> handleHit(Structure structure, String hit) {
        try {
            String chainId = hit.split("-")[0].split("_")[1];
            Chain chain = structure.select()
                    .chainName(chainId)
                    .asChain();
            String[] split = hit.split(",");
            double rmsd = Double.valueOf(split[1]);
            double pvalue = Double.valueOf(split[2]);

            String resIdSection = split[0].substring(5);
            List<AminoAcid> aminoAcids = Pattern.compile("_").splitAsStream(resIdSection)
                    .map(id -> id.split("-")[1])
                    .map(IdentifierFactory::createResidueIdentifier)
                    .map(residueIdentifier -> chain.select()
                            .residueIdentifier(residueIdentifier)
                            .asAminoAcid())
                    .collect(Collectors.toList());

            return Optional.of(hit.split(",")[0] + "," +
                    StandardFormat.format(rmsd) + "," +
                    StandardFormat.format(pvalue) + "," +
                    aminoAcids.stream()
                            .map(A01_AnnotateFeatures::composeAminoAcidFeatures)
                            .collect(Collectors.joining(","))
                    + ",?");
        } catch (Exception e) {
            e.printStackTrace();
            return Optional.empty();
        }
    }

    private static String composeAminoAcidFeatures(AminoAcid aminoAcid) {
        return StandardFormat.format(aminoAcid.getFeature(EnergyProfile.class).getSolvationEnergy()) + "," +
                StandardFormat.format(aminoAcid.getFeature(AccessibleSurfaceArea.class).getRelativeAccessibleSurfaceArea()) + "," +
                StandardFormat.format(aminoAcid.getFeature(LoopFraction.class).getLoopFraction()) + "," +
                aminoAcid.getFeature(PLIPInteractionContainer.class).getHydrogenBonds().size() + "," +
                aminoAcid.getFeature(PLIPInteractionContainer.class).getHydrophobicInteractions().size();
    }
}