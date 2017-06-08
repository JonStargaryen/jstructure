package studies.membrane.pdbtm;

import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.motif.SequenceMotifContainer;
import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import de.bioforscher.jstructure.feature.sse.SecondaryStructure;
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
    public static Collector<AminoAcid, ?, OccurrenceSummary> toOccurrenceSummary() {
        return Collector.of(OccurrenceSummary::new,
                OccurrenceSummary::accept,
                OccurrenceSummary::combine,
                Function.identity(),
                Collector.Characteristics.CONCURRENT);
    }

    static class OccurrenceSummary implements Consumer<AminoAcid> {
        private int count = 0;
        //                                tm, ntm
        private int[] region = new int[] { 0, 0 };
        private int[] motifs = new int[SequenceMotifDefinition.values().length];
        //                                                     H, S, c
        private int[] secondaryStructureElements = new int[] { 0, 0, 0 };
        //                                           a, y, m, c, s, l, w
        private int[] interactionTypes = new int[] { 0, 0, 0, 0, 0, 0, 0 };
        //                                          b, s, m
        private int[] interatingParts = new int[] { 0, 0, 0 };
        private boolean normalized = false;

        OccurrenceSummary() {
            this.count = 0;
        }

        @Override
        public void accept(AminoAcid aminoAcid) {
            Chain chain = aminoAcid.getParentChain();
            Protein protein = chain.getParentProtein();
            SecondaryStructure secondaryStructure = aminoAcid.getFeatureContainer().getFeature(SecondaryStructure.class);
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
                            interatingParts[0]++;
                        } else if(plipInteraction.isSideChainInteraction()) {
                            interatingParts[1]++;
                        } else {
                            interatingParts[2]++;
                        }
                    });
        }

        OccurrenceSummary combine(OccurrenceSummary other) {
            count += other.count;
            region = addVectors(region, other.region);
            motifs = addVectors(motifs, other.motifs);
            secondaryStructureElements = addVectors(secondaryStructureElements, other.secondaryStructureElements);
            interactionTypes = addVectors(interactionTypes, other.interactionTypes);
            interatingParts = addVectors(interatingParts, other.interatingParts);
            return this;
        }

        private int[] addVectors(int[] vector1, int[] vector2) {
            int[] result = new int[vector1.length];
            for(int index = 0; index <= result.length; index++) {
                result[index] = vector1[index] + vector2[index];
            }
            return result;
        }

        @Override
        public String toString() {
            return getHeaderLine() + System.lineSeparator() + getOccurrenceLine();
        }

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
                    IntStream.of(interatingParts)
                            .mapToObj(String::valueOf)
                            .collect(Collectors.joining("\t"));
        }

        /**
         * Normalize occurrence numbers
         */
        private void normalize() {
            //TODO impl
        }

        public static String getHeaderLine() {
            return "aminoAcids\ttm\tntm\t" +
                    Stream.of(SequenceMotifDefinition.values())
                            .map(SequenceMotifDefinition::name)
                            .collect(Collectors.joining("\t")) + "\thelix\tsheet\tcoil\thalogen\thydrogen\t" +
                    "metal\tpiCation\tpiStacking\tsalt\twater\tbackbone\tsideChain\tmixed";
        }
    }
}
