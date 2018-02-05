package de.bioforscher.start2fold.collect;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceArea;
import de.bioforscher.jstructure.feature.energyprofile.EgorAgreement;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfile;
import de.bioforscher.jstructure.feature.graphs.ResidueTopologicPropertiesContainer;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.loopfraction.LoopFraction;
import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.aminoacid.Proline;
import de.bioforscher.start2fold.Start2FoldConstants;
import de.bioforscher.start2fold.model.Start2FoldResidueAnnotation;
import de.bioforscher.start2fold.model.vector.FeatureVector;
import de.bioforscher.start2fold.model.vector.RawFeatureVector;
import de.bioforscher.start2fold.model.vector.SmoothedFeatureVector;
import de.bioforscher.start2fold.parser.Start2FoldXmlParser;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.function.ToDoubleFunction;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

public class A05_WriteEarlyFoldingClassificationArff {
    private static final Logger logger = LoggerFactory.getLogger(A05_WriteEarlyFoldingClassificationArff.class);

    public static void main(String[] args) throws IOException {
        String output = Files.lines(Start2FoldConstants.PANCSA_LIST)
                .map(A05_WriteEarlyFoldingClassificationArff::handleLine)
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.joining(System.lineSeparator(),
                        "@RELATION fc" + System.lineSeparator() +
                                /*"@ATTRIBUTE pdb string" + System.lineSeparator() +
                                "@ATTRIBUTE chain string" + System.lineSeparator() +
                                "@ATTRIBUTE res string" + System.lineSeparator() +
                                "@ATTRIBUTE aa string" + System.lineSeparator() +
                                "@ATTRIBUTE sse string" + System.lineSeparator() +*/
                                "@ATTRIBUTE sseSize numeric" + System.lineSeparator() +

                                "@ATTRIBUTE plip_l_hydrogen numeric" + System.lineSeparator() +
                                "@ATTRIBUTE plip_l_hydrophobic numeric" + System.lineSeparator() +
                                "@ATTRIBUTE plip_l_bb numeric" + System.lineSeparator() +
                                "@ATTRIBUTE plip_l_total numeric" + System.lineSeparator() +

                                "@ATTRIBUTE plip_nl_hydrogen numeric" + System.lineSeparator() +
                                "@ATTRIBUTE plip_nl_hydrophobic numeric" + System.lineSeparator() +
                                "@ATTRIBUTE plip_nl_bb numeric" + System.lineSeparator() +
                                "@ATTRIBUTE plip_nl_total numeric" + System.lineSeparator() +

                                "@ATTRIBUTE energy numeric" + System.lineSeparator() +
                                "@ATTRIBUTE egor numeric" + System.lineSeparator() +

//                                "@ATTRIBUTE eccount numeric" + System.lineSeparator() +
//                                "@ATTRIBUTE cumstrength numeric" + System.lineSeparator() +
//                                "@ATTRIBUTE ecstrength numeric" + System.lineSeparator() +
//                                "@ATTRIBUTE conservation numeric" + System.lineSeparator() +

                                "@ATTRIBUTE asa numeric" + System.lineSeparator() +
                                "@ATTRIBUTE loopFraction numeric" + System.lineSeparator() +

                                "@ATTRIBUTE plip_betweenness numeric" + System.lineSeparator() +
                                "@ATTRIBUTE plip_closeness numeric" + System.lineSeparator() +
                                "@ATTRIBUTE plip_clusteringcoefficient numeric" + System.lineSeparator() +
                                "@ATTRIBUTE plip_hydrogen_betweenness numeric" + System.lineSeparator() +
                                "@ATTRIBUTE plip_hydrogen_closeness numeric" + System.lineSeparator() +
                                "@ATTRIBUTE plip_hydrogen_clusteringcoefficient numeric" + System.lineSeparator() +
                                "@ATTRIBUTE plip_hydrophobic_betweenness numeric" + System.lineSeparator() +
                                "@ATTRIBUTE plip_hydrophobic_closeness numeric" + System.lineSeparator() +
                                "@ATTRIBUTE plip_hydrophobic_clusteringcoefficient numeric" + System.lineSeparator() +
                                "@ATTRIBUTE conv_betweenness numeric" + System.lineSeparator() +
                                "@ATTRIBUTE conv_closeness numeric" + System.lineSeparator() +
                                "@ATTRIBUTE conv_clusteringcoefficient numeric" + System.lineSeparator() +

                                "@ATTRIBUTE plip_distinct_neighborhoods numeric" + System.lineSeparator() +
                                "@ATTRIBUTE conv_distinct_neighborhoods numeric" + System.lineSeparator() +

                                "@ATTRIBUTE class {early,late}" + System.lineSeparator() +
                                "@DATA" + System.lineSeparator(),
                        ""));

        Start2FoldConstants.write(Start2FoldConstants.STATISTICS_DIRECTORY.resolve("foldingcores.arff"),
                output);
    }

    private static Optional<String> handleLine(String line) {
        try {
            logger.info("handling {}",
                    line);
            String[] split = line.split(";");
            String entryId = split[0];
            String pdbId = split[1];
            List<Integer> experimentIds = Pattern.compile(",")
                    .splitAsStream(split[2].replaceAll("\\[", "").replaceAll("]", ""))
                    .map(Integer::valueOf)
                    .collect(Collectors.toList());

            Structure structure = StructureParser.source(pdbId).parse();
            Chain chain = structure.getFirstChain();

            Start2FoldXmlParser.parseSpecificExperiment(chain,
                    Start2FoldConstants.XML_DIRECTORY.resolve(entryId + ".xml"),
                    experimentIds);

//            EvolutionaryCouplingParser.parseHotSpotFile(chain,
//                    Start2FoldConstants.COUPLING_DIRECTORY.resolve(entryId.toUpperCase() + "_hs.html"));

            List<AminoAcid> earlyFoldingResidues = chain.aminoAcids()
                    .filter(aminoAcid -> aminoAcid.getFeature(Start2FoldResidueAnnotation.class).isEarly())
                    .collect(Collectors.toList());

            List<AminoAcid> aminoAcids = chain.aminoAcids()
                    .collect(Collectors.toList());

            aminoAcids.forEach(aminoAcid -> {
                GenericSecondaryStructure sse = aminoAcid.getFeature(GenericSecondaryStructure.class);

                PLIPInteractionContainer plipInteractionContainer = aminoAcid.getFeature(PLIPInteractionContainer.class);
                PLIPInteractionContainer nonLocalPlipInteractionContainer = new PLIPInteractionContainer(null,
                        plipInteractionContainer
                                .getInteractions()
                                .stream()
                                // interactions have to be non-local
                                .filter(inter -> Math.abs(inter.getPartner1().getResidueIndex() - inter.getPartner2().getResidueIndex()) > 5)
                                .collect(Collectors.toList()));
                PLIPInteractionContainer localPlipInteractionContainer = new PLIPInteractionContainer(null,
                        plipInteractionContainer
                                .getInteractions()
                                .stream()
                                // interactions have to be local
                                .filter(inter -> !nonLocalPlipInteractionContainer.getInteractions().contains(inter))
                                .collect(Collectors.toList()));

                ResidueTopologicPropertiesContainer residueTopologicPropertiesContainer =
                        aminoAcid.getFeature(ResidueTopologicPropertiesContainer.class);

                // assign features to smooth
                RawFeatureVector featureVector = new RawFeatureVector(sse.getSurroundingSecondaryStructureElement(aminoAcid).getSize(),
                        localPlipInteractionContainer.getHydrogenBonds().size(),
                        localPlipInteractionContainer.getHydrophobicInteractions().size(),
                        localPlipInteractionContainer.getBackboneInteractions().size(),
                        localPlipInteractionContainer.getInteractions().size(),

                        nonLocalPlipInteractionContainer.getHydrogenBonds().size(),
                        nonLocalPlipInteractionContainer.getHydrophobicInteractions().size(),
                        nonLocalPlipInteractionContainer.getBackboneInteractions().size(),
                        nonLocalPlipInteractionContainer.getInteractions().size(),

                        aminoAcid.getFeature(EnergyProfile.class).getSolvationEnergy(),
                        aminoAcid.getFeature(EgorAgreement.class).getEgorPrediction(),

//                        aminoAcid.getFeature(HotSpotScoring.class).getEcCount(),
//                        aminoAcid.getFeature(HotSpotScoring.class).getCumStrength(),
//                        aminoAcid.getFeature(HotSpotScoring.class).getEcStrength(),
//                        aminoAcid.getFeature(HotSpotScoring.class).getConservation(),

                        aminoAcid.getFeature(AccessibleSurfaceArea.class).getRelativeAccessibleSurfaceArea(),

                        residueTopologicPropertiesContainer.getFullPlip().getBetweenness(),
                        residueTopologicPropertiesContainer.getFullPlip().getCloseness(),
                        residueTopologicPropertiesContainer.getFullPlip().getClusteringCoefficient(),
                        residueTopologicPropertiesContainer.getHydrogenPlip().getBetweenness(),
                        residueTopologicPropertiesContainer.getHydrogenPlip().getCloseness(),
                        residueTopologicPropertiesContainer.getHydrogenPlip().getClusteringCoefficient(),
                        residueTopologicPropertiesContainer.getHydrophobicPlip().getBetweenness(),
                        residueTopologicPropertiesContainer.getHydrophobicPlip().getCloseness(),
                        residueTopologicPropertiesContainer.getHydrophobicPlip().getClusteringCoefficient(),
                        residueTopologicPropertiesContainer.getConventional().getBetweenness(),
                        residueTopologicPropertiesContainer.getConventional().getCloseness(),
                        residueTopologicPropertiesContainer.getConventional().getClusteringCoefficient(),
                        residueTopologicPropertiesContainer.getFullPlip().getDistinctNeighborhoodCount(),
                       residueTopologicPropertiesContainer.getConventional().getDistinctNeighborhoodCount());

                aminoAcid.getFeatureContainer().addFeature(featureVector);
            });

            // smooth features
            aminoAcids.forEach(aminoAcid -> {
                int index = aminoAcid.getResidueIndex();
                List<AminoAcid> aminoAcidsToSmooth = new ArrayList<>();
                aminoAcidsToSmooth.add(aminoAcid);

                // get N-terminal residues
                for(int i = 0; i < 4; i++) {
                    int indexToGet = index - (i + 1);
                    if(indexToGet > -1) {
                        aminoAcidsToSmooth.add(aminoAcids.get(indexToGet));
                    }
                }

                // get C-terminal residues
                for(int i = 0; i < 4; i++) {
                    int indexToGet = index + (i + 1);
                    if(indexToGet < aminoAcids.size()) {
                        aminoAcidsToSmooth.add(aminoAcids.get(indexToGet));
                    }
                }

                // assign smoothed values
                smoothValues(aminoAcidsToSmooth, aminoAcid);
            });
            

            return Optional.of(aminoAcids.stream()
                    .filter(aminoAcid -> !(aminoAcid instanceof Proline))
                    .map(aminoAcid -> {
                        GenericSecondaryStructure sse = aminoAcid.getFeature(GenericSecondaryStructure.class);

                        SmoothedFeatureVector smoothedFeatureVector = aminoAcid.getFeature(SmoothedFeatureVector.class);

                        return /*pdbId + "," +
                                "A" + "," +
                                aminoAcid.getResidueIdentifier() + "," +
                                aminoAcid.getOneLetterCode() + "," +
                                sse.getSecondaryStructure().getOneLetterRepresentation() + "," +*/
                                StandardFormat.format(smoothedFeatureVector.getSecondaryStructureElementSize()) + "," +

                                StandardFormat.format(smoothedFeatureVector.getLocalHydrogen()) + "," +
                                StandardFormat.format(smoothedFeatureVector.getLocalHydrophobic()) + "," +
                                StandardFormat.format(smoothedFeatureVector.getLocalBackbone()) + "," +
                                StandardFormat.format(smoothedFeatureVector.getLocalInteractions()) + "," +

                                StandardFormat.format(smoothedFeatureVector.getNonLocalHydrogen()) + "," +
                                StandardFormat.format(smoothedFeatureVector.getNonLocalHydrophobic()) + "," +
                                StandardFormat.format(smoothedFeatureVector.getNonLocalBackbone()) + "," +
                                StandardFormat.format(smoothedFeatureVector.getNonLocalInteractions()) + "," +

                                StandardFormat.format(smoothedFeatureVector.getEnergy()) + "," +
                                StandardFormat.format(smoothedFeatureVector.getEgor()) + "," +

//                                StandardFormat.format(smoothedFeatureVector.getEccount()) + "," +
//                                StandardFormat.format(smoothedFeatureVector.getCumstrength()) + "," +
//                                StandardFormat.format(smoothedFeatureVector.getEcstrength()) + "," +
//                                StandardFormat.format(smoothedFeatureVector.getConservation()) + "," +

                                StandardFormat.format(smoothedFeatureVector.getRasa()) + "," +
                                StandardFormat.format(aminoAcid.getFeature(LoopFraction.class).getLoopFraction()) + "," + // already smoothed

                                StandardFormat.format(smoothedFeatureVector.getBetweenness()) + "," +
                                StandardFormat.format(smoothedFeatureVector.getCloseness()) + "," +
                                StandardFormat.format(smoothedFeatureVector.getClusteringCoefficient()) + "," +
                                StandardFormat.format(smoothedFeatureVector.getHydrogenBetweenness()) + "," +
                                StandardFormat.format(smoothedFeatureVector.getHydrogenCloseness()) + "," +
                                StandardFormat.format(smoothedFeatureVector.getHydrogenClusteringCoefficient()) + "," +
                                StandardFormat.format(smoothedFeatureVector.getHydrophobicBetweenness()) + "," +
                                StandardFormat.format(smoothedFeatureVector.getHydrophobicCloseness()) + "," +
                                StandardFormat.format(smoothedFeatureVector.getHydrophobicClusteringCoefficient()) + "," +
                                StandardFormat.format(smoothedFeatureVector.getConvBetweenness()) + "," +
                                StandardFormat.format(smoothedFeatureVector.getConvCloseness()) + "," +
                                StandardFormat.format(smoothedFeatureVector.getConvClusteringCoefficient()) + "," +

                                StandardFormat.format(smoothedFeatureVector.getDistinctNeighborhoods()) + "," +
                                StandardFormat.format(smoothedFeatureVector.getConvDistinctNeighborhoods()) + "," +

                                (earlyFoldingResidues.contains(aminoAcid) ? "early" : "late");
                    })
                    .collect(Collectors.joining(System.lineSeparator())));
        } catch (Exception e) {
            logger.warn("computation for {} failed",
                    line,
                    e);
            return Optional.empty();
        }
    }

    private static void smoothValues(List<AminoAcid> aminoAcidsToSmooth, AminoAcid aminoAcid) {
        double secondaryStructureElementSize = smoothValue(aminoAcidsToSmooth, FeatureVector::getSecondaryStructureElementSize);

        double localHydrogen = smoothValue(aminoAcidsToSmooth, FeatureVector::getLocalHydrogen);
        double localHydrophobic = smoothValue(aminoAcidsToSmooth, FeatureVector::getLocalHydrophobic);
        double localBackbone = smoothValue(aminoAcidsToSmooth, FeatureVector::getLocalBackbone);
        double localInteractions = smoothValue(aminoAcidsToSmooth, FeatureVector::getLocalInteractions);

        double nonLocalHydrogen = smoothValue(aminoAcidsToSmooth, FeatureVector::getNonLocalHydrogen);
        double nonLocalHydrophobic = smoothValue(aminoAcidsToSmooth, FeatureVector::getNonLocalHydrophobic);
        double nonLocalBackbone = smoothValue(aminoAcidsToSmooth, FeatureVector::getNonLocalBackbone);
        double nonLocalInteractions = smoothValue(aminoAcidsToSmooth, FeatureVector::getNonLocalInteractions);

        double energy = smoothValue(aminoAcidsToSmooth, FeatureVector::getEnergy);
        double egor = smoothValue(aminoAcidsToSmooth, FeatureVector::getEgor);

//        double eccount = smoothValue(aminoAcidsToSmooth, FeatureVector::getEccount);
//        double cumstrength = smoothValue(aminoAcidsToSmooth, FeatureVector::getCumstrength);
//        double ecstrength = smoothValue(aminoAcidsToSmooth, FeatureVector::getEcstrength);
//        double conservation = smoothValue(aminoAcidsToSmooth, FeatureVector::getConservation);

        double rasa = smoothValue(aminoAcidsToSmooth, FeatureVector::getRasa);

        double betweenness = smoothValue(aminoAcidsToSmooth, FeatureVector::getBetweenness);
        double closeness = smoothValue(aminoAcidsToSmooth, FeatureVector::getCloseness);
        double clusteringCoefficient = smoothValue(aminoAcidsToSmooth, FeatureVector::getClusteringCoefficient);

        double hydrogenBetweenness = smoothValue(aminoAcidsToSmooth, FeatureVector::getHydrogenBetweenness);
        double hydrogenCloseness = smoothValue(aminoAcidsToSmooth, FeatureVector::getHydrogenCloseness);
        double hydrogenClusteringCoefficient = smoothValue(aminoAcidsToSmooth, FeatureVector::getHydrogenClusteringCoefficient);

        double hydrophobicBetweenness = smoothValue(aminoAcidsToSmooth, FeatureVector::getHydrophobicBetweenness);
        double hydrophobicCloseness = smoothValue(aminoAcidsToSmooth, FeatureVector::getHydrophobicCloseness);
        double hydrophobicClusteringCoefficient = smoothValue(aminoAcidsToSmooth, FeatureVector::getHydrophobicClusteringCoefficient);

        double convBetweenness = smoothValue(aminoAcidsToSmooth, FeatureVector::getConvBetweenness);
        double convCloseness = smoothValue(aminoAcidsToSmooth, FeatureVector::getConvCloseness);
        double convClusteringCoefficient = smoothValue(aminoAcidsToSmooth, FeatureVector::getConvClusteringCoefficient);

        double distinctNeighborhoods = smoothValue(aminoAcidsToSmooth, FeatureVector::getDistinctNeighborhoods);
        double convDistinctNeighborhoods = smoothValue(aminoAcidsToSmooth, FeatureVector::getConvDistinctNeighborhoods);

        SmoothedFeatureVector smoothedFeatureVector = new SmoothedFeatureVector(secondaryStructureElementSize,

                localHydrogen,
                localHydrophobic,
                localBackbone,
                localInteractions,

                nonLocalHydrogen,
                nonLocalHydrophobic,
                nonLocalBackbone,
                nonLocalInteractions,

                energy,
                egor,

//                eccount,
//                cumstrength,
//                ecstrength,
//                conservation,

                rasa,

                betweenness,
                closeness,
                clusteringCoefficient,

                hydrogenBetweenness,
                hydrogenCloseness,
                hydrogenClusteringCoefficient,

                hydrophobicBetweenness,
                hydrophobicCloseness,
                hydrophobicClusteringCoefficient,

                convBetweenness,
                convCloseness,
                convClusteringCoefficient,

                distinctNeighborhoods,
                convDistinctNeighborhoods);

        aminoAcid.getFeatureContainer().addFeature(smoothedFeatureVector);
    }

    private static double smoothValue(List<AminoAcid> aminoAcidsToSmooth, ToDoubleFunction<FeatureVector> mapping) {
        return aminoAcidsToSmooth.stream()
                .map(aminoAcid -> aminoAcid.getFeature(RawFeatureVector.class))
                .mapToDouble(mapping)
                .average()
                .orElseThrow(() -> new IllegalArgumentException("could not compute average"));
    }
}
