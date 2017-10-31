package de.bioforscher.jstructure.membrane.contact;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionType;
import de.bioforscher.jstructure.feature.interactions.PLIPIntraMolecularAnnotator;
import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.feature.topology.OrientationsOfProteinsInMembranesAnnotator;
import de.bioforscher.jstructure.feature.topology.Topology;
import de.bioforscher.jstructure.membrane.MembraneConstants;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.GroupPrototype;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.jsoup.Jsoup;

import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class A03_AminoAcidPropensityForHelixHelixInteraction {
    private static final OrientationsOfProteinsInMembranesAnnotator ORIENTATIONS_OF_PROTEINS_IN_MEMBRANES_ANNOTATOR =
            new OrientationsOfProteinsInMembranesAnnotator();
    private static final PLIPIntraMolecularAnnotator PLIP_INTRA_MOLECULAR_ANNOTATOR =
            new PLIPIntraMolecularAnnotator();

    public static void main(String[] args) {
        String header = "aa;type;topology;count;rel" + System.lineSeparator();

        Path directory = MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY;

        Map<GroupPrototype, Map<String, Integer>> transmembraneCount = AminoAcid.Family.canonicalAminoAcids()
                .collect(Collectors.toMap(AminoAcid.Family::getGroupPrototype, family -> new HashMap<String, Integer>()));
        Map<GroupPrototype, Map<String, Integer>> nonTransmembraneCount = AminoAcid.Family.canonicalAminoAcids()
                .collect(Collectors.toMap(AminoAcid.Family::getGroupPrototype, family -> new HashMap<String, Integer>()));

        MembraneConstants.lines(directory.resolve("ids.list"))
                .forEach(line -> {
                    try {
                        System.out.println(line);
                        String pdbId = line.split("_")[0];
                        String chainId = line.split("_")[1];

                        Structure structure = StructureParser.source(directory.resolve("pdb").resolve(pdbId + ".pdb"))
                                .minimalParsing(true)
                                .parse();
                        Chain chain = structure.select()
                                .chainName(chainId)
                                .asChain();

                        ORIENTATIONS_OF_PROTEINS_IN_MEMBRANES_ANNOTATOR.process(structure,
                                Jsoup.parse(MembraneConstants.lines(directory.resolve("opm").resolve(pdbId + ".opm"))
                                        .collect(Collectors.joining(System.lineSeparator()))));
                        PLIP_INTRA_MOLECULAR_ANNOTATOR.process(chain,
                                Jsoup.parse(MembraneConstants.lines(directory.resolve("plip").resolve(line + ".plip"))
                                        .collect(Collectors.joining(System.lineSeparator()))));

                        chain.aminoAcids().forEach(aminoAcid -> {
                            boolean tm = aminoAcid.getFeature(Topology.class).isTransmembrane();
                            GenericSecondaryStructure.SecondaryStructureElement surroundingSecondaryStructureElement =
                                    aminoAcid.getFeature(GenericSecondaryStructure.class).getSurroundingSecondaryStructureElement(aminoAcid);
                            // skip non-transmembrane parts when there is no helix or it is too short
                            if (!tm && (!surroundingSecondaryStructureElement.getReducedType().equals("H") ||
                                    surroundingSecondaryStructureElement.getSize() <= MembraneConstants.MINIMUM_HELIX_LENGTH)) {
                                return;
                            }

                            Map<GroupPrototype, Map<String, Integer>> map = tm ? transmembraneCount : nonTransmembraneCount;

                            // count amino acid occurrence
                            Map<String, Integer> aminoAcidMap = map.get(aminoAcid.getGroupPrototype());
                            aminoAcidMap.merge("TOTAL", 1, Integer::sum);

                            // evaluate interaction types
                            PLIPInteractionContainer interactionContainer = aminoAcid.getFeature(PLIPInteractionContainer.class);
                            if(tm) {
                                interactionContainer.getInteractions()
                                        .stream()
                                        .filter(MembraneConstants::isTmHelixInteraction)
                                        .forEach(interaction -> aminoAcidMap.merge(interaction.getClass().getSimpleName(), 1, Integer::sum));
                            } else {
                                interactionContainer.getInteractions()
                                        .stream()
                                        .filter(MembraneConstants::isNtmHelixInteraction)
                                        .forEach(interaction -> aminoAcidMap.merge(interaction.getClass().getSimpleName(), 1, Integer::sum));
                            }
                        });
                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                });
        MembraneConstants.write(MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY.resolve("contact")
                        .resolve("aminoacids-interactions-relative-occurrence.csv"),
                header +
                        composeOutput(transmembraneCount, "I") + System.lineSeparator() +
                        composeOutput(nonTransmembraneCount, "o"));
    }

    private static String composeOutput(Map<GroupPrototype, Map<String, Integer>> map, String topology) {
        return AminoAcid.Family.canonicalAminoAcids()
                .map(family -> {
                    Map<String, Integer> aminoAcidMap = map.get(family.getGroupPrototype());
                    return Stream.of(PLIPInteractionType.values())
                            .map(plipInteractionType -> {
                                String type = plipInteractionType.getDescribingClass().getSimpleName();
                                int count = aminoAcidMap.getOrDefault(type, 0);
                                return family.getOneLetterCode() + ";" +
                                        type + ";" +
                                        topology + ";" +
                                        count + ";" +
                                        StandardFormat.format((count / (double) aminoAcidMap.getOrDefault("TOTAL", 1)) * 100);
                            })
                            .collect(Collectors.joining(System.lineSeparator()));
                })
                .collect(Collectors.joining(System.lineSeparator()));
    }
}
