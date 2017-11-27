package de.bioforscher.equant.train.collect;

import de.bioforscher.equant.EquantConstants;
import de.bioforscher.equant.train.TargetFunctionParser;
import de.bioforscher.jstructure.feature.energyprofile.EgorAgreementCalculator;
import de.bioforscher.jstructure.feature.evolution.EvolutionaryInformation;
import de.bioforscher.jstructure.feature.interactions.PLIPIntraMolecularAnnotator;
import de.bioforscher.jstructure.feature.loopfraction.LoopFractionCalculator;
import de.bioforscher.jstructure.feature.rigidity.DynaMineBridge;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.jsoup.Jsoup;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.nio.file.Path;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

public class T06_WriteTrainingArff {
    private static final Logger logger = LoggerFactory.getLogger(T06_WriteTrainingArff.class);
    private static final Path DIRECTORY = EquantConstants.CASP9_DIRECTORY;
    private static final Path DYNAMINE_DIRECTORY = DIRECTORY.resolve("dynamine");
    private static final Path PLIP_DIRECTORY = DIRECTORY.resolve("plip");
    private static final Path PREDICTIONS_DIRECTORY = DIRECTORY.resolve("predictions");
    private static final Path PSIBLAST_DIRECTORY = DIRECTORY.resolve("psiblast");
    private static final Path RESULTS_DIRECTORY = DIRECTORY.resolve("results");
    private static final TargetFunctionParser TARGET_FUNCTION_PARSER = new TargetFunctionParser();
    private static final DynaMineBridge DYNA_MINE_BRIDGE = new DynaMineBridge();
    private static final PLIPIntraMolecularAnnotator PLIP_INTRA_MOLECULAR_ANNOTATOR = new PLIPIntraMolecularAnnotator();
    private static final LoopFractionCalculator LOOP_FRACTION_CALCULATOR = new LoopFractionCalculator();
    private static final EgorAgreementCalculator EGOR_AGREEMENT_CALCULATOR = new EgorAgreementCalculator();

    public static void main(String[] args) {
        String output = EquantConstants.list(PREDICTIONS_DIRECTORY)
                .map(T06_WriteTrainingArff::handleTarget)
                .collect(Collectors.joining(System.lineSeparator()));

        System.out.println(output);
    }

    private static String handleTarget(Path path) {
        String targetName = path.toFile().getName();

        logger.info("handling target {}", targetName);

        List<EvolutionaryInformation> psiBlastProfile = EquantConstants.lines(PSIBLAST_DIRECTORY.resolve(targetName + ".psiblast"))
                .map(Double::valueOf)
                .map(value -> new EvolutionaryInformation(null, null, value))
                .collect(Collectors.toList());

        return EquantConstants.list(path)
                .map(model -> handleModel(model, targetName))
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.joining(System.lineSeparator()));
    }

    private static Optional<String> handleModel(Path path,
                                                String targetName) {
        try {
            String modelName = path.toFile().getName();
            logger.info("handling {}", modelName);
            Structure model = StructureParser.source(path).parse();
            Chain chain = model.chains().findFirst().get();

            // parse target function scores
            TARGET_FUNCTION_PARSER.process(model, RESULTS_DIRECTORY.resolve(targetName).resolve(modelName + ".lga"));

            // parse dynamine scores
            DYNA_MINE_BRIDGE.process(chain,
                    EquantConstants.lines(DYNAMINE_DIRECTORY.resolve(targetName + ".dynamine"))
                            .collect(Collectors.joining(System.lineSeparator())));

            // parse PLIP data
            PLIP_INTRA_MOLECULAR_ANNOTATOR.process(chain,
                    Jsoup.parse(EquantConstants.lines(PLIP_DIRECTORY.resolve(targetName).resolve(modelName + ".plip"))
                            .collect(Collectors.joining(System.lineSeparator()))));

            // compute egor agreement
            EGOR_AGREEMENT_CALCULATOR.process(model);

            // compute loop fraction scores
            LOOP_FRACTION_CALCULATOR.process(model);

            return Optional.of(model.aminoAcids()
                    .map(aminoAcid -> handleResidue(aminoAcid, modelName, targetName))
                    .filter(Optional::isPresent)
                    .map(Optional::get)
                    .collect(Collectors.joining(System.lineSeparator())));
        } catch (Exception e) {
            e.printStackTrace();
            return Optional.empty();
        }
    }

    private static Optional<String> handleResidue(AminoAcid aminoAcid,
                                                  String modelName,
                                                  String targetName) {
        try {
            TargetFunctionParser.TargetFunction targetFunction = aminoAcid.getFeature(TargetFunctionParser.TargetFunction.class);

            return Optional.of(modelName + "," +
                    targetName + "," +
                    aminoAcid.getResidueIdentifier() + "," +
                    targetFunction.getDistance() + "," +
                    targetFunction.getSscore() + ",");
        } catch (Exception e) {
            e.printStackTrace();
            return Optional.empty();
        }
    }
}
