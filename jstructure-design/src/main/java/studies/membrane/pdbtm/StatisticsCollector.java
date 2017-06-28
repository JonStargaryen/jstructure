package studies.membrane.pdbtm;

import de.bioforscher.jstructure.feature.interactions.*;
import de.bioforscher.jstructure.feature.motif.SequenceMotifContainer;
import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import de.bioforscher.jstructure.feature.sse.dssp.DSSPSecondaryStructure;
import de.bioforscher.jstructure.feature.topology.MembraneContainer;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.util.function.Consumer;
import java.util.function.Function;
import java.util.stream.Collector;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Custom statistic collectors.
 * Created by bittrich on 6/7/17.
 */
public class StatisticsCollector {
    public static Collector<PLIPInteraction, ?, InteractingAminoAcidSummary> toInteractingAminoAcidSummary() {
        return Collector.of(InteractingAminoAcidSummary::new,
                InteractingAminoAcidSummary::accept,
                InteractingAminoAcidSummary::combine,
                Function.identity(),
                Collector.Characteristics.CONCURRENT);
    }

    interface CustomCollector {
        String getOccurrenceLine();

        String getHeaderLine();

        default String getLine() {
            return getHeaderLine() + System.lineSeparator() + getOccurrenceLine();
        }

        default int[] addVectors(int[] vector1, int[] vector2) {
            int[] result = new int[vector1.length];
            for(int index = 0; index <= result.length; index++) {
                result[index] = vector1[index] + vector2[index];
            }
            return result;
        }
    }

    public static class InteractingAminoAcidSummary implements Consumer<PLIPInteraction>, CustomCollector {
        private int count = 0;
        private int[] motifs = new int[SequenceMotifDefinition.values().length];
        //                                           a, y, m, c, s, l, w
        private int[] interactionTypes = new int[] { 0, 0, 0, 0, 0, 0, 0 };
        //                                          b, s, m
        private int[] interactingParts = new int[] { 0, 0, 0 };
        private boolean normalized = false;

        @Override
        public void accept(PLIPInteraction plipInteraction) {
            SequenceMotifContainer sequenceMotifContainer = plipInteraction
                    .getPartner1()
                    .getParentChain()
                    .getParentProtein()
                    .getFeatureContainer()
                    .getFeature(SequenceMotifContainer.class);

            count++;

            // sequence motifs
            sequenceMotifContainer.getEmbeddingSequenceMotifsFor(plipInteraction.getPartner1()).getSequenceMotifs()
                    .forEach(sequenceMotif -> motifs[sequenceMotif.getMotifDefinition().ordinal()]++);
            sequenceMotifContainer.getEmbeddingSequenceMotifsFor(plipInteraction.getPartner2()).getSequenceMotifs()
                    .forEach(sequenceMotif -> motifs[sequenceMotif.getMotifDefinition().ordinal()]++);

            // interaction type
            if(plipInteraction instanceof HalogenBond) {
                interactionTypes[0]++;
            } else if(plipInteraction instanceof HydrogenBond) {
                interactionTypes[1]++;
            } else if(plipInteraction instanceof MetalComplex) {
                interactionTypes[2]++;
            } else if(plipInteraction instanceof PiCationInteraction) {
                interactionTypes[3]++;
            } else if(plipInteraction instanceof PiStacking) {
                interactionTypes[4]++;
            } else if(plipInteraction instanceof SaltBridge) {
                interactionTypes[5]++;
            } else {
                interactionTypes[6]++;
            }

            // interacting atoms
            if(plipInteraction.isBackboneInteraction()) {
                interactingParts[0]++;
            } else if(plipInteraction.isSideChainInteraction()) {
                interactingParts[1]++;
            } else {
                interactingParts[2]++;
            }
        }

        InteractingAminoAcidSummary combine(InteractingAminoAcidSummary other) {
            count += other.count;
            motifs = addVectors(motifs, other.motifs);
            interactionTypes = addVectors(interactionTypes, other.interactionTypes);
            interactingParts = addVectors(interactingParts, other.interactingParts);
            return this;
        }

        @Override
        public String getOccurrenceLine() {
            if(!normalized) {
                normalize();
            }

            return count + "\t" +
                    IntStream.of(motifs)
                            .mapToObj(String::valueOf)
                            .collect(Collectors.joining("\t")) + "\t" +
                    IntStream.of(interactionTypes)
                            .mapToObj(String::valueOf)
                            .collect(Collectors.joining("\t")) + "\t" +
                    IntStream.of(interactingParts)
                            .mapToObj(String::valueOf)
                            .collect(Collectors.joining("\t"));
        }

        /**
         * Normalize occurrence numbers, e.g. by the number of possible interactions or sequence motifs
         */
        private void normalize() {
            //TODO impl
        }

        @Override
        public String getHeaderLine() {
            return "aminoAcids\t" +
                    Stream.of(SequenceMotifDefinition.values())
                            .map(SequenceMotifDefinition::name)
                            .collect(Collectors.joining("\t")) + "\thalogen\thydrogen\t" +
                    "metal\tpiCation\tpiStacking\tsalt\twater\tbackbone\tsideChain\tmixed";
        }
    }

    public static Collector<AminoAcid, ?, AminoAcidSummary> toAminoAcidSummary() {
        return Collector.of(AminoAcidSummary::new,
                AminoAcidSummary::accept,
                AminoAcidSummary::combine,
                Function.identity(),
                Collector.Characteristics.CONCURRENT);
    }

    public static class AminoAcidSummary implements Consumer<AminoAcid>, CustomCollector {
        private int count = 0;
        //                                tm, ntm
        private int[] region = new int[] { 0, 0 };
        private int[] motifs = new int[SequenceMotifDefinition.values().length];
        //                                                     H, S, c
        private int[] secondaryStructureElements = new int[] { 0, 0, 0 };
        //                                           a, y, m, c, s, l, w
        private int[] interactionTypes = new int[] { 0, 0, 0, 0, 0, 0, 0 };
        //                                          b, s, m
        private int[] interactingParts = new int[] { 0, 0, 0 };
        private boolean normalized = false;

        public AminoAcidSummary() {
            this.count = 0;
        }

        @Override
        public void accept(AminoAcid aminoAcid) {
            Chain chain = aminoAcid.getParentChain();
            Protein protein = chain.getParentProtein();
            DSSPSecondaryStructure secondaryStructure = aminoAcid.getFeatureContainer().getFeature(DSSPSecondaryStructure.class);
            PLIPInteractionContainer plipInteractionContainer = chain.getFeatureContainer().getFeature(PLIPInteractionContainer.class);
            MembraneContainer membraneContainer = protein.getFeatureContainer().getFeature(MembraneContainer.class);
            SequenceMotifContainer sequenceMotifContainer = protein.getFeatureContainer().getFeature(SequenceMotifContainer.class);

            count++;

            // region
            if(membraneContainer.isTransmembraneGroup(aminoAcid)) {
                region[0]++;
            } else {
                region[1]++;
            }

            // sequence motifs
            sequenceMotifContainer.getEmbeddingSequenceMotifsFor(aminoAcid).getSequenceMotifs()
                    .forEach(sequenceMotif -> motifs[sequenceMotif.getMotifDefinition().ordinal()]++);

            // secondary structure
            switch(secondaryStructure.getSecondaryStructure().getReducedRepresentation()) {
                case "H":
                    secondaryStructureElements[0]++;
                    break;
                case "E":
                    secondaryStructureElements[1]++;
                    break;
                default:
                    secondaryStructureElements[2]++;
            }

            PLIPInteractionContainer groupSpecificInteractionContainer = plipInteractionContainer.getInteractionsFor(aminoAcid);
            // interactions
            interactionTypes[0] += groupSpecificInteractionContainer.getHalogenBonds().size();
            interactionTypes[1] += groupSpecificInteractionContainer.getHydrogenBonds().size();
            interactionTypes[2] += groupSpecificInteractionContainer.getMetalComplexes().size();
            interactionTypes[3] += groupSpecificInteractionContainer.getPiCationInteractions().size();
            interactionTypes[4] += groupSpecificInteractionContainer.getPiStackings().size();
            interactionTypes[5] += groupSpecificInteractionContainer.getSaltBridges().size();
            interactionTypes[6] += groupSpecificInteractionContainer.getWaterBridges().size();

            // interacting atoms
            groupSpecificInteractionContainer.getInteractions()
                    .forEach(plipInteraction -> {
                        if(plipInteraction.isBackboneInteraction()) {
                            interactingParts[0]++;
                        } else if(plipInteraction.isSideChainInteraction()) {
                            interactingParts[1]++;
                        } else {
                            interactingParts[2]++;
                        }
                    });
        }

        AminoAcidSummary combine(AminoAcidSummary other) {
            count += other.count;
            region = addVectors(region, other.region);
            motifs = addVectors(motifs, other.motifs);
            secondaryStructureElements = addVectors(secondaryStructureElements, other.secondaryStructureElements);
            interactionTypes = addVectors(interactionTypes, other.interactionTypes);
            interactingParts = addVectors(interactingParts, other.interactingParts);
            return this;
        }

        @Override
        public String getOccurrenceLine() {
            if(!normalized) {
                normalize();
            }

            return count + "\t" +
                    IntStream.of(region)
                            .mapToObj(String::valueOf)
                            .collect(Collectors.joining("\t")) + "\t" +
                    IntStream.of(motifs)
                            .mapToObj(String::valueOf)
                            .collect(Collectors.joining("\t")) + "\t" +
                    IntStream.of(secondaryStructureElements)
                            .mapToObj(String::valueOf)
                            .collect(Collectors.joining("\t")) + "\t" +
                    IntStream.of(interactionTypes)
                            .mapToObj(String::valueOf)
                            .collect(Collectors.joining("\t")) + "\t" +
                    IntStream.of(interactingParts)
                            .mapToObj(String::valueOf)
                            .collect(Collectors.joining("\t"));
        }

        /**
         * Normalize occurrence numbers, e.g. by the number of possible interactions or sequence motifs
         */
        private void normalize() {
            //TODO impl
        }

        @Override
        public String getHeaderLine() {
            return "aminoAcids\ttm\tntm\t" +
                    Stream.of(SequenceMotifDefinition.values())
                            .map(SequenceMotifDefinition::name)
                            .collect(Collectors.joining("\t")) + "\thelix\tsheet\tcoil\thalogen\thydrogen\t" +
                    "metal\tpiCation\tpiStacking\tsalt\twater\tbackbone\tsideChain\tmixed";
        }
    }
}
