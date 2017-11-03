package de.bioforscher.jstructure.membrane;

import de.bioforscher.jstructure.feature.interactions.PLIPInteraction;
import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.feature.topology.MembraneContainer;
import de.bioforscher.jstructure.feature.topology.Topology;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.selection.IntegerRange;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.UncheckedIOException;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Optional;
import java.util.stream.Stream;

/**
 * Shared class of project constants and convenience functions used in the membrane modules.
 */
public class MembraneConstants {
    private static final Logger logger = LoggerFactory.getLogger(MembraneConstants.class);
    public static final Path GIT_DIRECTORY = Paths.get(System.getProperty("user.home")).resolve("git");
    public static final Path PHD_DIRECTORY = GIT_DIRECTORY.resolve("phd_sb_repo");
    public static final Path DATA_DIRECTORY = PHD_DIRECTORY.resolve("data");
    public static final Path DATASETS_DIRECTORY = DATA_DIRECTORY.resolve("datasets");
    // test performance and optimal setup of module annotation
    public static final Path MODULARITY_DATASET_DIRECTORY = DATASETS_DIRECTORY.resolve("modularity");
    public static final Path PDBTM_NR_ALPHA_DATASET_DIRECTORY = DATASETS_DIRECTORY.resolve("pdbtm_nr_alpha");
    public static final Path PDBTM_NR_ALPHA_DATASET_PDB_DIRECTORY = PDBTM_NR_ALPHA_DATASET_DIRECTORY.resolve("pdb");
    public static final Path PDBTM_NR_ALPHA_DATASET_OPM_DIRECTORY = PDBTM_NR_ALPHA_DATASET_DIRECTORY.resolve("opm");
    public static final Path PDBTM_NR_ALPHA_DATASET_PLIP_DIRECTORY = PDBTM_NR_ALPHA_DATASET_DIRECTORY.resolve("plip");
    public static final Path PDBTM_NR_ALPHA_DATASET_NETWORK_DIRECTORY = PDBTM_NR_ALPHA_DATASET_DIRECTORY.resolve("network");
    public static final Path PDBTM_NR_BETA_DATASET_DIRECTORY = DATASETS_DIRECTORY.resolve("pdbtm_nr_beta");
    public static final Path PDBTM_NR_BETA_DATASET_PDB_DIRECTORY = PDBTM_NR_BETA_DATASET_DIRECTORY.resolve("pdb");
    public static final Path PDBTM_NR_BETA_DATASET_OPM_DIRECTORY = PDBTM_NR_BETA_DATASET_DIRECTORY.resolve("opm");
    public static final Path PDBTM_NR_BETA_DATASET_PLIP_DIRECTORY = PDBTM_NR_BETA_DATASET_DIRECTORY.resolve("plip");
    public static final Path PDBTM_NR_BETA_DATASET_NETWORK_DIRECTORY = PDBTM_NR_BETA_DATASET_DIRECTORY.resolve("network");
    public static final Path FOLDING_CORES_DIRECTORY = DATASETS_DIRECTORY.resolve("foldingcores");
    public static final Path DIVISION_DIRECTORY = DATASETS_DIRECTORY.resolve("division");
    public static final int MINIMUM_HELIX_LENGTH = 15;

    public static Stream<String> lines(Path path) {
        try {
            return Files.lines(path);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    public static Stream<String> lines(String url) {
        try {
            return new BufferedReader(new InputStreamReader(new URL(url).openStream())).lines();
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    public static Stream<Path> list(Path path) {
        try {
            return Files.list(path);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    public static void write(Path path, String content) {
        ensureDirectoriesExist(path);
        try {
            Files.write(path, content.getBytes());
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    private static void ensureDirectoriesExist(Path path) {
        if(Files.isDirectory(path)) {
            throw new IllegalArgumentException(path.toFile().getAbsolutePath() + " is no regular file - cannot process");
        }

        Path parentDirectory = path.getParent();
        if(!Files.exists(parentDirectory)) {
            logger.info("creating directory '{}'", parentDirectory.toFile().getAbsolutePath());
            try {
                Files.createDirectories(parentDirectory);
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
        }
    }

    public static void move(Path inPath, Path outPath) {
        try {
            Files.move(inPath, outPath);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    public static boolean isNtmHelixInteraction(PLIPInteraction plipInteraction) {
        if(plipInteraction.getPartner1().getFeature(Topology.class).isTransmembrane() ||
                plipInteraction.getPartner2().getFeature(Topology.class).isTransmembrane()) {
            return false;
        }

        AminoAcid aminoAcid1 = (AminoAcid) plipInteraction.getPartner1();
        AminoAcid aminoAcid2 = (AminoAcid) plipInteraction.getPartner2();


        GenericSecondaryStructure.SecondaryStructureElement secondaryStructureElement1 =
                aminoAcid1.getFeature(GenericSecondaryStructure.class).getSurroundingSecondaryStructureElement(aminoAcid1);
        GenericSecondaryStructure.SecondaryStructureElement secondaryStructureElement2 =
                aminoAcid2.getFeature(GenericSecondaryStructure.class).getSurroundingSecondaryStructureElement(aminoAcid2);

        return secondaryStructureElement1.getReducedType().equals("H") &&
                secondaryStructureElement2.getReducedType().equals("H") &&
                secondaryStructureElement1.getSize() > MembraneConstants.MINIMUM_HELIX_LENGTH &&
                secondaryStructureElement2.getSize() > MembraneConstants.MINIMUM_HELIX_LENGTH &&
                !secondaryStructureElement1.equals(secondaryStructureElement2);

    }

    public static boolean isTmHelixInteraction(PLIPInteraction plipInteraction) {
        if(!plipInteraction.getPartner1().getFeature(Topology.class).isTransmembrane() ||
                !plipInteraction.getPartner2().getFeature(Topology.class).isTransmembrane()) {
            return false;
        }

        Structure structure = plipInteraction.getPartner1().getParentStructure();
        MembraneContainer membraneContainer = structure.getFeature(MembraneContainer.class);
        Optional<IntegerRange> segment1 = membraneContainer.getEmbeddingTransmembraneSegment(plipInteraction.getPartner1());
        Optional<IntegerRange> segment2 = membraneContainer.getEmbeddingTransmembraneSegment(plipInteraction.getPartner2());

        return segment1.isPresent() && segment2.isPresent() && !segment1.equals(segment2);
    }

    private static final int SEQUENCE_WINDOW = 10;
    public static String surroundingSequence(AminoAcid aminoAcid) {
        StringBuilder upstreamSequence = new StringBuilder();
        StringBuilder downstreamSequence = new StringBuilder();

        Optional<AminoAcid> aminoAcidOptional = aminoAcid.getPreviousAminoAcid();
        for(int i = 0; i < SEQUENCE_WINDOW; i++) {
            if(aminoAcidOptional.isPresent()) {
                upstreamSequence.append(aminoAcidOptional.get().getOneLetterCode());
                aminoAcidOptional = aminoAcidOptional.get().getPreviousAminoAcid();
            } else {
                upstreamSequence.append("-");
                aminoAcidOptional = Optional.empty();
            }
        }

        aminoAcidOptional = aminoAcid.getNextAminoAcid();
        for(int i = 0; i < SEQUENCE_WINDOW; i++) {
            if(aminoAcidOptional.isPresent()) {
                downstreamSequence.append(aminoAcidOptional.get().getOneLetterCode());
                aminoAcidOptional = aminoAcidOptional.get().getNextAminoAcid();
            } else {
                downstreamSequence.append("-");
                aminoAcidOptional = Optional.empty();
            }
        }

        return upstreamSequence.reverse().toString() + aminoAcid.getOneLetterCode() + downstreamSequence;
    }

    public static double minMaxNormalize(double v, double min, double max) {
        return (v - min) / (max  - min);
    }
}
