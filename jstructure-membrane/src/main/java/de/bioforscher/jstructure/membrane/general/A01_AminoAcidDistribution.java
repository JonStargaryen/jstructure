package de.bioforscher.jstructure.membrane.general;

import de.bioforscher.jstructure.StandardFormat;
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

/**
 * Frequencies of individual amino acids in the membrane resp. in helices of length >15 outside the membrane.
 */
public class A01_AminoAcidDistribution {
    private static final OrientationsOfProteinsInMembranesAnnotator orientationsOfProteinsInMembranesAnnotator =
            new OrientationsOfProteinsInMembranesAnnotator();

    public static void main(String[] args) {
        String sse = "H";
        //TODO literature - length must be greater than 15
        int minimalHelixLength = MembraneConstants.MINIMUM_HELIX_LENGTH;
        Path directory = MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY;

        Map<GroupPrototype, Integer> transmembraneCount = new HashMap<>();
        Map<GroupPrototype, Integer> nonTransmembraneCount = new HashMap<>();

        MembraneConstants.lines(directory.resolve("ids.list"))
                .forEach(id -> {
                    System.out.println(id);
                    String pdbId = id.split("_")[0];
                    String chainId = id.split("_")[1];

                    try {
                        Structure structure = StructureParser.source(directory.resolve("pdb").resolve(pdbId + ".pdb"))
                                .minimalParsing(true)
                                .parse();
                        Chain chain = structure.select()
                                .chainId(chainId)
                                .asChain();
                        orientationsOfProteinsInMembranesAnnotator.process(structure,
                                Jsoup.parse(MembraneConstants.lines(directory.resolve("opm").resolve(pdbId + ".opm"))
                                        .collect(Collectors.joining(System.lineSeparator()))));

                        chain.aminoAcids().forEach(aminoAcid -> {
                            boolean tm = aminoAcid.getFeature(Topology.class).isTransmembrane();
                            GenericSecondaryStructure.SecondaryStructureElement surroundingSecondaryStructureElement =
                                    aminoAcid.getFeature(GenericSecondaryStructure.class).getSurroundingSecondaryStructureElement(aminoAcid);
                            // skip non-transmembrane parts when there are no helix or too short
                            if (!tm && (!surroundingSecondaryStructureElement.getReducedType().equals(sse) ||
                                    surroundingSecondaryStructureElement.getSize() <= minimalHelixLength)) {
                                return;
                            }

                            Map<GroupPrototype, Integer> map = tm ? transmembraneCount : nonTransmembraneCount;
                            // initializes with 1 if key is absent, otherwise increment by 1
                            map.merge(aminoAcid.getGroupPrototype(), 1, Integer::sum);
                        });
                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                });

        double tmTotal = transmembraneCount.values()
                .stream()
                .mapToDouble(Double::valueOf)
                .sum();
        double ntmTotal = nonTransmembraneCount.values()
                .stream()
                .mapToDouble(Double::valueOf)
                .sum();

        MembraneConstants.write(directory.resolve("general").resolve("aminoacids-frequencies.csv"),
                "aa;topology;count;rel" + System.lineSeparator() +
                        AminoAcid.Family.canonicalAminoAcids()
                                .map(family -> family.getOneLetterCode() + ";" +
                                        "I;" +
                                        transmembraneCount.get(family.getGroupPrototype()) + ";" +
                                        StandardFormat.format(transmembraneCount.get(family.getGroupPrototype()) / tmTotal * 100))
                                .collect(Collectors.joining(System.lineSeparator())) + System.lineSeparator() +
                        AminoAcid.Family.canonicalAminoAcids()
                                .map(family -> family.getOneLetterCode() + ";" +
                                        "o;" +
                                        nonTransmembraneCount.get(family.getGroupPrototype()) + ";" +
                                        StandardFormat.format(nonTransmembraneCount.get(family.getGroupPrototype()) / ntmTotal * 100))
                                .collect(Collectors.joining(System.lineSeparator())));
    }
}
