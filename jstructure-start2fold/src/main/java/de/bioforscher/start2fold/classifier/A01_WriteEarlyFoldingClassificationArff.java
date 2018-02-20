package de.bioforscher.start2fold.classifier;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.feature.loopfraction.LoopFraction;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.aminoacid.Proline;
import de.bioforscher.start2fold.Start2FoldConstants;
import de.bioforscher.start2fold.model.Start2FoldResidueAnnotation;
import de.bioforscher.start2fold.model.vector.RawFeatureVector;
import de.bioforscher.start2fold.model.vector.SmoothedFeatureVector;
import de.bioforscher.start2fold.parser.Start2FoldXmlParser;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Optional;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

public class A01_WriteEarlyFoldingClassificationArff {
    private static final Logger logger = LoggerFactory.getLogger(A01_WriteEarlyFoldingClassificationArff.class);

    public static void main(String[] args) throws IOException {
        String output = Files.lines(Start2FoldConstants.PANCSA_LIST)
                .map(A01_WriteEarlyFoldingClassificationArff::handleLine)
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.joining(System.lineSeparator(),
                        "@RELATION fc" + System.lineSeparator() +

                                "@ATTRIBUTE energy numeric" + System.lineSeparator() +
                                "@ATTRIBUTE egor numeric" + System.lineSeparator() +

                                "@ATTRIBUTE sse_size numeric" + System.lineSeparator() +
                                "@ATTRIBUTE loop_fraction numeric" + System.lineSeparator() +

                                "@ATTRIBUTE rasa numeric" + System.lineSeparator() +

                                "@ATTRIBUTE plip_local_contacts numeric" + System.lineSeparator() +
                                "@ATTRIBUTE plip_local_hbonds numeric" + System.lineSeparator() +
                                "@ATTRIBUTE plip_local_hydrophobic numeric" + System.lineSeparator() +
                                "@ATTRIBUTE plip_local_backbone numeric" + System.lineSeparator() +

                                "@ATTRIBUTE plip_long_range_contacts numeric" + System.lineSeparator() +
                                "@ATTRIBUTE plip_long_range_hbonds numeric" + System.lineSeparator() +
                                "@ATTRIBUTE plip_long_range_hydrophobic numeric" + System.lineSeparator() +
                                "@ATTRIBUTE plip_long_range_backbone numeric" + System.lineSeparator() +

                                "@ATTRIBUTE plip_betweenness numeric" + System.lineSeparator() +
                                "@ATTRIBUTE plip_closeness numeric" + System.lineSeparator() +
                                "@ATTRIBUTE plip_clusteringcoefficient numeric" + System.lineSeparator() +

                                "@ATTRIBUTE plip_hbonds_betweenness numeric" + System.lineSeparator() +
                                "@ATTRIBUTE plip_hbonds_closeness numeric" + System.lineSeparator() +
                                "@ATTRIBUTE plip_hbonds_clusteringcoefficient numeric" + System.lineSeparator() +

                                "@ATTRIBUTE plip_hydrophobic_betweenness numeric" + System.lineSeparator() +
                                "@ATTRIBUTE plip_hydrophobic_closeness numeric" + System.lineSeparator() +
                                "@ATTRIBUTE plip_hydrophobic_clusteringcoefficient numeric" + System.lineSeparator() +

                                "@ATTRIBUTE conv_betweenness numeric" + System.lineSeparator() +
                                "@ATTRIBUTE conv_closeness numeric" + System.lineSeparator() +
                                "@ATTRIBUTE conv_clusteringcoefficient numeric" + System.lineSeparator() +

                                "@ATTRIBUTE plip_neighborhoods numeric" + System.lineSeparator() +
                                "@ATTRIBUTE conv_neighborhoods numeric" + System.lineSeparator() +

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

            Structure structure = StructureParser.fromPdbId(pdbId).parse();
            Chain chain = structure.getFirstChain();

            Start2FoldXmlParser.parseSpecificExperiment(chain,
                    Start2FoldConstants.XML_DIRECTORY.resolve(entryId + ".xml"),
                    experimentIds);

            List<AminoAcid> earlyFoldingResidues = chain.aminoAcids()
                    .filter(aminoAcid -> aminoAcid.getFeature(Start2FoldResidueAnnotation.class).isEarly())
                    .collect(Collectors.toList());

            List<AminoAcid> aminoAcids = chain.aminoAcids()
                    .collect(Collectors.toList());

            aminoAcids.forEach(RawFeatureVector::assignRawFeatureVector);

            // smooth features
            aminoAcids.forEach(aminoAcid -> SmoothedFeatureVector.assignSmoothedFeatureVector(aminoAcids, aminoAcid));
            

            return Optional.of(aminoAcids.stream()
                    .filter(aminoAcid -> !(aminoAcid instanceof Proline))
                    .map(aminoAcid -> {
                        SmoothedFeatureVector smoothedFeatureVector = aminoAcid.getFeature(SmoothedFeatureVector.class);

                        return StandardFormat.format(smoothedFeatureVector.getEnergy()) + "," +
                                StandardFormat.format(smoothedFeatureVector.getEgor()) + "," +

                                StandardFormat.format(smoothedFeatureVector.getSecondaryStructureElementSize()) + "," +
                                StandardFormat.format(aminoAcid.getFeature(LoopFraction.class).getLoopFraction()) + "," + // already smoothed

                                StandardFormat.format(smoothedFeatureVector.getRasa()) + "," +

                                StandardFormat.format(smoothedFeatureVector.getLocalInteractions()) + "," +
                                StandardFormat.format(smoothedFeatureVector.getLocalHydrogen()) + "," +
                                StandardFormat.format(smoothedFeatureVector.getLocalHydrophobic()) + "," +
                                StandardFormat.format(smoothedFeatureVector.getLocalBackbone()) + "," +

                                StandardFormat.format(smoothedFeatureVector.getNonLocalInteractions()) + "," +
                                StandardFormat.format(smoothedFeatureVector.getNonLocalHydrogen()) + "," +
                                StandardFormat.format(smoothedFeatureVector.getNonLocalHydrophobic()) + "," +
                                StandardFormat.format(smoothedFeatureVector.getNonLocalBackbone()) + "," +

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
}
