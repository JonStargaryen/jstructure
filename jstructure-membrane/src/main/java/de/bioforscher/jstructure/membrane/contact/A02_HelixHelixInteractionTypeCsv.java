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
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import org.jsoup.Jsoup;

import java.nio.file.Path;
import java.util.stream.Collectors;

public class A02_HelixHelixInteractionTypeCsv {
    private static final OrientationsOfProteinsInMembranesAnnotator ORIENTATIONS_OF_PROTEINS_IN_MEMBRANES_ANNOTATOR =
            new OrientationsOfProteinsInMembranesAnnotator();
    private static final PLIPIntraMolecularAnnotator PLIP_INTRA_MOLECULAR_ANNOTATOR =
            new PLIPIntraMolecularAnnotator();

    public static void main(String[] args) {
        String header = "type;topology;count;rel" + System.lineSeparator();

        Path directory = MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY;
        int[] tm = new int[9];
        int[] ntm = new int[9];
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
                            PLIPInteractionContainer interactionContainer = aminoAcid.getFeature(PLIPInteractionContainer.class);
                            if(aminoAcid.getFeature(Topology.class).isTransmembrane()) {
                                tm[8]++;
                                interactionContainer.getInteractions()
                                        .stream()
                                        .filter(MembraneConstants::isTmHelixInteraction)
                                        .forEach(interaction -> {
                                            tm[PLIPInteractionType.valueOf(interaction.getClass()).get().ordinal()]++;
                                        });
                            } else if(aminoAcid.getFeature(GenericSecondaryStructure.class).getSecondaryStructure().isHelixType()) {
                                ntm[8]++;
                                interactionContainer.getInteractions()
                                        .stream()
                                        .filter(MembraneConstants::isNtmHelixInteraction)
                                        .forEach(interaction -> {
                                            ntm[PLIPInteractionType.valueOf(interaction.getClass()).get().ordinal()]++;
                                        });
                            }
                        });
                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                });
        MembraneConstants.write(MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY.resolve("contact").resolve("interactions-relative-occurrence.csv"),
                header +
                PLIPInteractionType.HALOGEN_BOND + ";I;" + tm[0] + ";" + StandardFormat.format((tm[0] / (double) tm[8]) * 100) + System.lineSeparator() +
                PLIPInteractionType.HYDROGEN_BOND + ";I;" + tm[1] + ";" + StandardFormat.format((tm[1] / (double) tm[8]) * 100) + System.lineSeparator() +
                PLIPInteractionType.HYDROPHOBIC_INTERACTION + ";I;" + tm[2] + ";" + StandardFormat.format((tm[2] / (double) tm[8]) * 100) + System.lineSeparator() +
                PLIPInteractionType.METAL_COMPLEX + ";I;" + tm[3] + ";" + StandardFormat.format((tm[3] / (double) tm[8]) * 100) + System.lineSeparator() +
                PLIPInteractionType.PI_CATION_INTERACTION + ";I;" + tm[4] + ";" + StandardFormat.format((tm[4] / (double) tm[8]) * 100) + System.lineSeparator() +
                PLIPInteractionType.PI_STACKING + ";I;" + tm[5] + ";" + StandardFormat.format((tm[5] / (double) tm[8]) * 100) + System.lineSeparator() +
                PLIPInteractionType.SALT_BRIDGE + ";I;" + tm[6] + ";" + StandardFormat.format((tm[6] / (double) tm[8]) * 100) + System.lineSeparator() +
                PLIPInteractionType.WATER_BRIDGE + ";I;" + tm[7] + ";" + StandardFormat.format((tm[7] / (double) tm[8]) * 100) + System.lineSeparator() +
                PLIPInteractionType.HALOGEN_BOND + ";o;" + ntm[0] + ";" + StandardFormat.format((ntm[0] / (double) ntm[8]) * 100) + System.lineSeparator() +
                PLIPInteractionType.HYDROGEN_BOND + ";o;" + ntm[1] + ";" + StandardFormat.format((ntm[1] / (double) ntm[8]) * 100) + System.lineSeparator() +
                PLIPInteractionType.HYDROPHOBIC_INTERACTION + ";o;" + ntm[2] + ";" + StandardFormat.format((ntm[2] / (double) ntm[8]) * 100) + System.lineSeparator() +
                PLIPInteractionType.METAL_COMPLEX + ";o;" + ntm[3] + ";" + StandardFormat.format((ntm[3] / (double) ntm[8]) * 100) + System.lineSeparator() +
                PLIPInteractionType.PI_CATION_INTERACTION + ";o;" + ntm[4] + ";" + StandardFormat.format((ntm[4] / (double) ntm[8]) * 100) + System.lineSeparator() +
                PLIPInteractionType.PI_STACKING + ";o;" + ntm[5] + ";" + StandardFormat.format((ntm[5] / (double) ntm[8]) * 100) + System.lineSeparator() +
                PLIPInteractionType.SALT_BRIDGE + ";o;" + ntm[6] + ";" + StandardFormat.format((ntm[6] / (double) ntm[8]) * 100) + System.lineSeparator() +
                PLIPInteractionType.WATER_BRIDGE + ";o;" + ntm[7] + ";" + StandardFormat.format((ntm[7] / (double) ntm[8]) * 100) + System.lineSeparator());
    }
}
