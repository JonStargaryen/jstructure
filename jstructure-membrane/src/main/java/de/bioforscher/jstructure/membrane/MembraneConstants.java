package de.bioforscher.jstructure.membrane;

import de.bioforscher.jstructure.feature.evolution.EvolutionaryInformation;
import de.bioforscher.jstructure.feature.interactions.PLIPInteraction;
import de.bioforscher.jstructure.feature.interactions.PLIPIntraMolecularAnnotator;
import de.bioforscher.jstructure.feature.rigidity.DynaMineBridge;
import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.feature.topology.MembraneContainer;
import de.bioforscher.jstructure.feature.topology.OrientationsOfProteinsInMembranesAnnotator;
import de.bioforscher.jstructure.feature.topology.Topology;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.selection.IntegerRange;
import de.bioforscher.testutil.FileUtils;
import org.jsoup.Jsoup;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

/**
 * Shared class of project constants and convenience functions used in the membrane modules.
 */
public class MembraneConstants extends FileUtils {
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
    public static final Path COUPLING_DIRECTORY = DATASETS_DIRECTORY.resolve("membrane_couplings");
    public static final int MINIMUM_HELIX_LENGTH = 15;

    private static final OrientationsOfProteinsInMembranesAnnotator ORIENTATIONS_OF_PROTEINS_IN_MEMBRANES_ANNOTATOR =
            new OrientationsOfProteinsInMembranesAnnotator();
    private static final PLIPIntraMolecularAnnotator PLIP_INTRA_MOLECULAR_ANNOTATOR =
            new PLIPIntraMolecularAnnotator();
    private static final DynaMineBridge DYNA_MINE_BRIDGE =
            new DynaMineBridge();

    public static Optional<WrappedChain> handleLine(String id, Path directory) {
        try {
            logger.info("annotating {}",
                    id);
            String pdbId = id.split("_")[0];
            String chainId = id.split("_")[1];

            Structure structure = StructureParser.source(directory.resolve("pdb").resolve(pdbId + ".pdb"))
                    .minimalParsing(true)
                    .parse();
            Chain chain = structure.select()
                    .chainName(chainId)
                    .asChain();

            ORIENTATIONS_OF_PROTEINS_IN_MEMBRANES_ANNOTATOR.process(structure,
                    Jsoup.parse(MembraneConstants.lines(directory.resolve("opm").resolve(pdbId + ".opm"))
                            .collect(Collectors.joining(System.lineSeparator()))));

            // skip chains with no or a single tm region
//            int tmRegions = 0;
//            boolean inTmRegion = false;
//            List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());
//            for (int i = 0; i < aminoAcids.size(); i++) {
//                boolean tm = aminoAcids.get(i).getFeature(Topology.class).isTransmembrane();
//                if(i == 0) {
//                    inTmRegion = tm;
//                    continue;
//                }
//
//                if(inTmRegion && !tm) {
//                    tmRegions++;
//                }
//
//                inTmRegion = tm;
//            }
//            if(inTmRegion) {
//                tmRegions++;
//            }
            int index = -1;
            List<List<AminoAcid>> transmembraneHelices = new ArrayList<>();
            boolean wasPreviouslyTransmembrane = false;
            List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());
            for(AminoAcid aminoAcid : aminoAcids) {
                boolean transmembrane = aminoAcid.getFeature(Topology.class).isTransmembrane();
                if(transmembrane) {
                    if(wasPreviouslyTransmembrane) {
                        transmembraneHelices.get(index).add(aminoAcid);
                    } else {
                        index++;
                        wasPreviouslyTransmembrane = true;
                        List<AminoAcid> region = new ArrayList<>();
                        region.add(aminoAcid);
                        transmembraneHelices.add(region);
                    }
                } else {
                    if(wasPreviouslyTransmembrane) {
                        transmembraneHelices.add(new ArrayList<>());
                    }
                    wasPreviouslyTransmembrane = false;
                }
            }
            if(wasPreviouslyTransmembrane) {
                transmembraneHelices.get(index).add(aminoAcids.get(aminoAcids.size() - 1));
            }
            transmembraneHelices.removeIf(Collection::isEmpty);

            if(transmembraneHelices.size() < 2) {
                logger.warn("skipping structure with {} TM-regions",
                        transmembraneHelices.size());
                return Optional.empty();
            }

            PLIP_INTRA_MOLECULAR_ANNOTATOR.process(chain,
                    Jsoup.parse(MembraneConstants.lines(directory.resolve("plip").resolve(id + ".plip"))
                            .collect(Collectors.joining(System.lineSeparator()))));
            DYNA_MINE_BRIDGE.process(chain,
                    MembraneConstants.lines(directory.resolve("dynamine").resolve(id + "_backbone.pred"))
                            .collect(Collectors.joining(System.lineSeparator())));
            List<Double> evolutionaryScores = MembraneConstants.lines(directory.resolve("evol").resolve(id + ".evol"))
                    .mapToDouble(Double::valueOf)
                    .boxed()
                    .collect(Collectors.toList());
            boolean evolScoresSane = aminoAcids.size() == evolutionaryScores.size();
            if(evolScoresSane) {
                for(int i = 0; i < aminoAcids.size(); i++) {
                    aminoAcids.get(i).getFeatureContainer().addFeature(new EvolutionaryInformation(null, null, evolutionaryScores.get(i)));
                }
            }

            // assign kink data
            boolean kinkException = false;
            try {
                MembraneConstants.lines(directory.resolve("kinks").resolve(pdbId + ".kinks"))
                        .filter(kinkLine -> !kinkLine.startsWith("pdb_code"))
                        .filter(kinkLine -> kinkLine.substring(4, 5).equals(chainId))
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
            } catch (Exception e) {
                e.printStackTrace();
                kinkException = true;
            }
            boolean kinkFinderDataSane = !kinkException;

            return Optional.of(new WrappedChain(chain,
                    pdbId,
                    chainId,
                    evolScoresSane,
                    kinkFinderDataSane,
                    transmembraneHelices));
        } catch (Exception e) {
            logger.warn("compution failed for {}",
                    id,
                    e);
            return Optional.empty();
        }
    }

    public static class WrappedChain {
        private final Chain chain;
        private final String pdbId;
        private final String chainId;
        private final boolean evolScoresSane;
        private final boolean kinkFinderDataSane;
        private final List<List<AminoAcid>> transmembraneHelices;

        WrappedChain(Chain chain,
                     String pdbId,
                     String chainId,
                     boolean evolScoresSane,
                     boolean kinkFinderDataSane,
                     List<List<AminoAcid>> transmembraneHelices) {
            this.chain = chain;
            this.pdbId = pdbId;
            this.chainId = chainId;
            this.evolScoresSane = evolScoresSane;
            this.kinkFinderDataSane = kinkFinderDataSane;
            this.transmembraneHelices = transmembraneHelices;
        }

        public Chain getChain() {
            return chain;
        }

        public String getPdbId() {
            return pdbId;
        }

        public String getChainId() {
            return chainId;
        }

        public boolean isEvolScoresSane() {
            return evolScoresSane;
        }

        public boolean isKinkFinderDataSane() {
            return kinkFinderDataSane;
        }

        public List<List<AminoAcid>> getTransmembraneHelices() {
            return transmembraneHelices;
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
}
