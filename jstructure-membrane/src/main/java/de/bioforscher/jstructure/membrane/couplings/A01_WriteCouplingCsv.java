package de.bioforscher.jstructure.membrane.couplings;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.feature.interactions.PLIPInteraction;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.interactions.PLIPIntraMolecularAnnotator;
import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.feature.topology.OrientationsOfProteinsInMembranesAnnotator;
import de.bioforscher.jstructure.feature.topology.Topology;
import de.bioforscher.jstructure.membrane.MembraneConstants;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.jsoup.Jsoup;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

/**
 * For a given list: annotate all interactions. Assess them by their evolutionary coupling score.
 */
public class A01_WriteCouplingCsv {
    private static final Logger logger = LoggerFactory.getLogger(A01_WriteCouplingCsv.class);
    private static final Path directory = MembraneConstants.COUPLING_DIRECTORY;
    private static final PLIPIntraMolecularAnnotator PLIP_INTRA_MOLECULAR_ANNOTATOR =
            new PLIPIntraMolecularAnnotator();
    private static final OrientationsOfProteinsInMembranesAnnotator ORIENTATIONS_OF_PROTEINS_IN_MEMBRANES_ANNOTATOR =
            new OrientationsOfProteinsInMembranesAnnotator();

    public static void main(String[] args) {
        String output = MembraneConstants.lines(directory.resolve("ids.list"))
                .map(A01_WriteCouplingCsv::handleLine)
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.joining(System.lineSeparator()));

        MembraneConstants.write(directory.resolve("results").resolve("couplings.csv"),
                output);
    }

    private static Optional<String> handleLine(String id) {
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

            String output = chain.getFeature(PLIPInteractionContainer.class)
                    .getInteractions()
                    .stream()
                    .filter(MembraneConstants::isTmHelixInteraction)
                    .map(interaction -> handleInteraction(pdbId,
                            chainId,
                            interaction))
                    .filter(Optional::isPresent)
                    .map(Optional::get)
                    .collect(Collectors.joining(System.lineSeparator()));

            if(output.isEmpty()) {
                return Optional.empty();
            }

            return Optional.of(output);
        } catch (Exception e) {
            logger.warn("compution failed for {}",
                    id,
                    e);
            return Optional.empty();
        }
    }

    private static Optional<String> handleInteraction(String pdbId,
                                                      String chainId,
                                                      PLIPInteraction interaction) {
        try {
            AminoAcid aminoAcid1 = (AminoAcid) interaction.getPartner1();
            AminoAcid aminoAcid2 = (AminoAcid) interaction.getPartner2();

            PLIPInteractionContainer interactions1 = aminoAcid1.getFeature(PLIPInteractionContainer.class);
            PLIPInteractionContainer interactions2 = aminoAcid2.getFeature(PLIPInteractionContainer.class);

            // evaluate interaction types
            PLIPInteractionContainer nonLocalInteractions1 = new PLIPInteractionContainer(null,
                    interactions1
                            .getInteractions()
                            .stream()
                            // interactions have to be non-local
                            .filter(inter -> Math.abs(inter.getPartner1().getResidueIdentifier().getResidueNumber() - inter.getPartner2().getResidueIdentifier().getResidueNumber()) > 6)
                            .collect(Collectors.toList()));
            PLIPInteractionContainer nonLocalInteractions2 = new PLIPInteractionContainer(null,
                    interactions2
                            .getInteractions()
                            .stream()
                            // interactions have to be non-local
                            .filter(inter -> Math.abs(inter.getPartner1().getResidueIdentifier().getResidueNumber() - inter.getPartner2().getResidueIdentifier().getResidueNumber()) > 6)
                            .collect(Collectors.toList()));

            Optional<Double> score = MembraneConstants.lines(MembraneConstants.COUPLING_DIRECTORY.resolve("couplings").resolve(pdbId + "_" + chainId + ".txt"))
                    .map(line -> line.split("\\s+"))
                    .filter(split -> (aminoAcid1.getResidueIdentifier().toString().equals(split[0]) && aminoAcid2.getResidueIdentifier().toString().equals(split[2])) ||
                            (aminoAcid1.getResidueIdentifier().toString().equals(split[2]) && aminoAcid2.getResidueIdentifier().toString().equals(split[0])))
                    .findFirst()
                    .map(split -> split[5])
                    .map(Double::valueOf);

            GenericSecondaryStructure secondaryStructure1 = aminoAcid1.getFeature(GenericSecondaryStructure.class);
            GenericSecondaryStructure secondaryStructure2 = aminoAcid2.getFeature(GenericSecondaryStructure.class);
            GenericSecondaryStructure.SecondaryStructureElement surroundingSecondaryStructureElement1 =
                    aminoAcid1.getFeature(GenericSecondaryStructure.class).getSurroundingSecondaryStructureElement(aminoAcid1);
            GenericSecondaryStructure.SecondaryStructureElement surroundingSecondaryStructureElement2 =
                    aminoAcid2.getFeature(GenericSecondaryStructure.class).getSurroundingSecondaryStructureElement(aminoAcid2);

            return Optional.of(pdbId + ";" +
                    chainId + ";" +
                    aminoAcid1.getResidueIdentifier() + ";" +
                    aminoAcid2.getResidueIdentifier() + ";" +
                    aminoAcid1.getOneLetterCode() + ";" +
                    aminoAcid2.getOneLetterCode() + ";" +
                    MembraneConstants.surroundingSequence(aminoAcid1) + ";" +
                    MembraneConstants.surroundingSequence(aminoAcid2) + ";" +
                    secondaryStructure1.getSecondaryStructure().getReducedRepresentation() + ";" +
                    secondaryStructure1.getSecondaryStructure().getOneLetterRepresentation() + ";" +
                    secondaryStructure2.getSecondaryStructure().getReducedRepresentation() + ";" +
                    secondaryStructure2.getSecondaryStructure().getOneLetterRepresentation() + ";" +
                    surroundingSecondaryStructureElement1.getSize() + ";" +
                    surroundingSecondaryStructureElement2.getSize() + ";" +

                    interaction.getClass().getSimpleName() + ";" +
                    StandardFormat.format(aminoAcid1.getCa().calculate().distance(aminoAcid2.getCa())) + ";" +

                    score.map(StandardFormat::format)
                            .orElse("NA") + ";" +

                    interactions1.getBackboneInteractions().size() + ";" +
                    interactions2.getBackboneInteractions().size() + ";" +

                    nonLocalInteractions1.getHalogenBonds().size() + ";" +
                    nonLocalInteractions1.getHydrogenBonds().size() + ";" +
                    nonLocalInteractions1.getHydrophobicInteractions().size() + ";" +
                    nonLocalInteractions1.getMetalComplexes().size() + ";" +
                    nonLocalInteractions1.getPiCationInteractions().size() + ";" +
                    nonLocalInteractions1.getPiStackings().size() + ";" +
                    nonLocalInteractions1.getSaltBridges().size() + ";" +
                    nonLocalInteractions1.getWaterBridges().size() + ";" +
                    nonLocalInteractions1.getInteractions().size() + ";" +

                    nonLocalInteractions2.getHalogenBonds().size() + ";" +
                    nonLocalInteractions2.getHydrogenBonds().size() + ";" +
                    nonLocalInteractions2.getHydrophobicInteractions().size() + ";" +
                    nonLocalInteractions2.getMetalComplexes().size() + ";" +
                    nonLocalInteractions2.getPiCationInteractions().size() + ";" +
                    nonLocalInteractions2.getPiStackings().size() + ";" +
                    nonLocalInteractions2.getSaltBridges().size() + ";" +
                    nonLocalInteractions2.getWaterBridges().size() + ";" +
                    nonLocalInteractions2.getInteractions().size());
        } catch (Exception e) {
            e.printStackTrace();
            return Optional.empty();
        }
    }
}
