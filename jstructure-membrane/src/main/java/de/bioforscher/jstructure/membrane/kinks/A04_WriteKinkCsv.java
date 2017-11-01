package de.bioforscher.jstructure.membrane.kinks;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.interactions.PLIPIntraMolecularAnnotator;
import de.bioforscher.jstructure.feature.rigidity.BackboneRigidity;
import de.bioforscher.jstructure.feature.rigidity.DynaMineBridge;
import de.bioforscher.jstructure.feature.topology.OrientationsOfProteinsInMembranesAnnotator;
import de.bioforscher.jstructure.feature.topology.Topology;
import de.bioforscher.jstructure.membrane.MembraneConstants;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.jsoup.Jsoup;

import java.nio.file.Path;
import java.util.Optional;
import java.util.stream.Collectors;

public class A04_WriteKinkCsv {
    private static final Path DIRECTORY = MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY;
    private static final OrientationsOfProteinsInMembranesAnnotator ORIENTATIONS_OF_PROTEINS_IN_MEMBRANES_ANNOTATOR =
            new OrientationsOfProteinsInMembranesAnnotator();
    private static final DynaMineBridge DYNA_MINE_BRIDGE = new DynaMineBridge();
    private static final PLIPIntraMolecularAnnotator PLIP_INTRA_MOLECULAR_ANNOTATOR = new PLIPIntraMolecularAnnotator();

    public static void main(String[] args) {
        String header = "pdbId;chainId;resId;aa;topology;contact;total;rigidity;kink;sign;angle" + System.lineSeparator();

        String output = MembraneConstants.lines(DIRECTORY.resolve("ids.list"))
                .map(A04_WriteKinkCsv::handleId)
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.joining(System.lineSeparator()));
        MembraneConstants.write(MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY.resolve("contact")
                        .resolve("aminoacids-rigidity-by-contact-and-kinks.csv"),
                header + output);
    }

    private static Optional<String> handleId(String id) {
        System.out.println(id);
        String pdbId = id.split("_")[0];
        String chainId = id.split("_")[1];

        try {
            Structure structure = StructureParser.source(DIRECTORY.resolve("pdb").resolve(pdbId + ".pdb"))
                    .minimalParsing(true)
                    .parse();
            Chain chain = structure.select()
                    .chainId(chainId)
                    .asChain();

            // assign kink data
            MembraneConstants.lines(DIRECTORY.resolve("kinks").resolve(pdbId + ".kinks"))
                    .filter(line -> !line.startsWith("pdb_code"))
                    .filter(line -> line.substring(4,5 ).equals(chainId))
                    .forEach(kinkLine -> {
                        String[] split = kinkLine.split(",");
                        int kinkPosition = Integer.valueOf(split[3]);
                        double kinkAngle = Double.valueOf(split[6]);
                        boolean significantKink = kinkAngle > 15;
                        chain.select()
                                .aminoAcids()
                                .residueNumber(kinkPosition)
                                .asAminoAcid()
                                .getFeatureContainer()
                                .addFeature(new Kink(kinkAngle, significantKink));
                    });

            ORIENTATIONS_OF_PROTEINS_IN_MEMBRANES_ANNOTATOR.process(structure,
                    Jsoup.parse(MembraneConstants.lines(DIRECTORY.resolve("opm").resolve(pdbId + ".opm"))
                            .collect(Collectors.joining(System.lineSeparator()))));
            PLIP_INTRA_MOLECULAR_ANNOTATOR.process(chain,
                    Jsoup.parse(MembraneConstants.lines(DIRECTORY.resolve("plip").resolve(id + ".plip"))
                            .collect(Collectors.joining(System.lineSeparator()))));
            DYNA_MINE_BRIDGE.process(chain,
                    MembraneConstants.lines(DIRECTORY.resolve("dynamine").resolve(id + "_backbone.pred"))
                            .collect(Collectors.joining(System.lineSeparator())));

            double min = chain.aminoAcids()
                    .mapToDouble(aa -> aa.getFeature(BackboneRigidity.class).getBackboneRigidity())
                    .min()
                    .getAsDouble();
            double max = chain.aminoAcids()
                    .mapToDouble(aa -> aa.getFeature(BackboneRigidity.class).getBackboneRigidity())
                    .max()
                    .getAsDouble();

            String output = chain.aminoAcids()
                    .map(aminoAcid -> handleAminoAcid(pdbId, chainId, aminoAcid, min, max))
                    .filter(Optional::isPresent)
                    .map(Optional::get)
                    .collect(Collectors.joining(System.lineSeparator()));

            return Optional.of(output);
        } catch (Exception e) {
            e.printStackTrace();
            return Optional.empty();
        }
    }

    static class Kink extends FeatureContainerEntry {
        private final double angle;
        private final boolean significantKink;

        public Kink(double angle, boolean significantKink) {
            super(null);
            this.angle = angle;
            this.significantKink = significantKink;
        }

        public double getAngle() {
            return angle;
        }

        public boolean isSignificantKink() {
            return significantKink;
        }
    }

    private static Optional<String> handleAminoAcid(String pdbId,
                                                    String chainId,
                                                    AminoAcid aminoAcid,
                                                    double min,
                                                    double max) {
        try {
            boolean tm = aminoAcid.getFeature(Topology.class).isTransmembrane();

            // evaluate interaction types
            PLIPInteractionContainer interactionContainer = aminoAcid.getFeature(PLIPInteractionContainer.class);
            boolean inContact = tm ?
                    interactionContainer.getInteractions()
                            .stream()
                            .anyMatch(MembraneConstants::isTmHelixInteraction) :
                    interactionContainer.getInteractions()
                            .stream()
                            .anyMatch(MembraneConstants::isNtmHelixInteraction);

            Optional<Kink> featureOptional = aminoAcid.getFeatureContainer().getFeatureOptional(Kink.class);
            boolean kinkPresent = featureOptional.isPresent();

            return Optional.of(pdbId + ";" +
                    chainId + ";" +
                    aminoAcid.getResidueIdentifier() + ";" +
                    aminoAcid.getOneLetterCode() + ";" +
                    (tm ? "I" : "o") + ";" +
                    (inContact ? "contact" : "no-contact") + ";" +
                    interactionContainer.getInteractions().size() + ";" +
                    StandardFormat.format(minMaxNormalize(aminoAcid.getFeature(BackboneRigidity.class).getBackboneRigidity(), min, max)) + ";" +
                    (kinkPresent ? "K" : "n") + ";" +
                    (kinkPresent ? (featureOptional.get().isSignificantKink() ? "S" : "i") : "NA") + ";" +
                    (kinkPresent ? featureOptional.get().getAngle() : "NA"));
        } catch (Exception e) {
            e.printStackTrace();
            return Optional.empty();
        }
    }

    private static double minMaxNormalize(double v, double min, double max) {
        return (v - min) / (max  - min);
    }
}
