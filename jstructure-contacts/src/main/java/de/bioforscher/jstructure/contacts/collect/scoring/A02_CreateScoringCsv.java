package de.bioforscher.jstructure.contacts.collect.scoring;

import de.bioforscher.jstructure.contacts.ContactsConstants;
import de.bioforscher.jstructure.feature.interactions.PLIPIntraMolecularAnnotator;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import org.jsoup.Jsoup;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Optional;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

public class A02_CreateScoringCsv {
    private static final Logger logger = LoggerFactory.getLogger(A02_CreateScoringCsv.class);
    private static final Path DIRECTORY = ContactsConstants.GIT_DIRECTORY.resolve("phd_sb_repo")
            .resolve("data")
            .resolve("equant")
            .resolve("casp")
            .resolve("CASP9")
            .resolve("predictions");
    private static final Path PLIP_DIRECTORY = DIRECTORY.getParent().resolve("plip");
    private static final PLIPIntraMolecularAnnotator PLIP_INTRA_MOLECULAR_ANNOTATOR = new PLIPIntraMolecularAnnotator();

    public static void main(String[] args) {
        String output = ContactsConstants.walk(DIRECTORY)
                .filter(Files::isRegularFile)
                .map(A02_CreateScoringCsv::handleModel)
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.joining(System.lineSeparator()));

        System.out.println(output);
    }

    private static Optional<String> handleModel(Path path) {
        String modelName = path.toFile().getName();
        String targetName = modelName.substring(0, 5);
        try {
            Path baseDirectory = path.getParent()
                    .getParent()
                    .getParent();
            Structure target = StructureParser.source(baseDirectory.resolve("targets")
                    .resolve(targetName + ".pdb"))
                    .parse();
            Structure model = StructureParser.source(path).parse();

            // parse PLIP data
            PLIP_INTRA_MOLECULAR_ANNOTATOR.process(target.chains().findFirst().get(),
                    Jsoup.parse(ContactsConstants.lines(PLIP_DIRECTORY.resolve(targetName).resolve(modelName + ".plip"))
                            .collect(Collectors.joining(System.lineSeparator()))));
            PLIP_INTRA_MOLECULAR_ANNOTATOR.process(model.chains().findFirst().get(),
                    Jsoup.parse(ContactsConstants.lines(PLIP_DIRECTORY.resolve(targetName).resolve(targetName + ".plip"))
                            .collect(Collectors.joining(System.lineSeparator()))));

            ContactMapScoring contactMapScoring = ContactMapScorer.score(target, model);

            LGAGlobalScores lgaGlobalScores = LGAGlobalParser.parseScores(baseDirectory.resolve("results")
                    .resolve(targetName)
                    .resolve(modelName + ".lga"));

            return Optional.of(targetName + "," +
                    modelName + "," +
                    lgaGlobalScores.getRmsd() + "," +
                    lgaGlobalScores.getGdt() + "," +
                    lgaGlobalScores.getLga_s3() + "," +
                    lgaGlobalScores.getLga_q() + "," +
                    contactMapScoring.getConventionalConfusionMatrix().getSensitivity() + "," +
                    contactMapScoring.getConventionalConfusionMatrix().getSpecificity() + "," +
                    contactMapScoring.getConventionalConfusionMatrix().getAccuracy() + "," +
                    contactMapScoring.getConventionalConfusionMatrix().getF1Score() + "," +
                    contactMapScoring.getConventionalConfusionMatrix().getMCC()
            );
        } catch (Exception e) {
            logger.warn("failed computation for {}",
                    modelName,
                    e);
            return Optional.empty();
        }
    }

    public static class LGAGlobalParser {
        public static LGAGlobalScores parseScores(Path path) {
            String resultLine = ContactsConstants.lines(path)
                    .filter(line -> line.startsWith("SUMMARY(GDT)"))
                    .findFirst()
                    .get();
            String[] split = Pattern.compile("\\s+").split(resultLine);
            return new LGAGlobalScores(
                    Double.valueOf(split[5]),
                    Double.valueOf(split[6]),
                    Double.valueOf(split[7]),
                    Double.valueOf(split[8])
            );
        }
    }

    public static class LGAGlobalScores {
        private final double rmsd;
        private final double gdt;
        private final double lga_s3;
        private final double lga_q;

        public LGAGlobalScores(double rmsd,
                               double gdt,
                               double lga_s3,
                               double lga_q) {
            this.rmsd = rmsd;
            this.gdt = gdt;
            this.lga_s3 = lga_s3;
            this.lga_q = lga_q;
        }

        public double getRmsd() {
            return rmsd;
        }

        public double getGdt() {
            return gdt;
        }

        public double getLga_s3() {
            return lga_s3;
        }

        public double getLga_q() {
            return lga_q;
        }
    }
}