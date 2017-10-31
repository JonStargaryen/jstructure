package de.bioforscher.jstructure.membrane.contact;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.feature.interactions.HydrogenBond;
import de.bioforscher.jstructure.feature.interactions.HydrophobicInteraction;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.interactions.PLIPIntraMolecularAnnotator;
import de.bioforscher.jstructure.feature.topology.OrientationsOfProteinsInMembranesAnnotator;
import de.bioforscher.jstructure.feature.topology.Topology;
import de.bioforscher.jstructure.mathematics.Pair;
import de.bioforscher.jstructure.membrane.MembraneConstants;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.jsoup.Jsoup;

import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;
import java.util.stream.Collectors;

public class A05_CooccurrenceMatrixCsv {
    private static final Path DIRECTORY = MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY;
    private static final OrientationsOfProteinsInMembranesAnnotator ORIENTATIONS_OF_PROTEINS_IN_MEMBRANES_ANNOTATOR =
            new OrientationsOfProteinsInMembranesAnnotator();
    private static final PLIPIntraMolecularAnnotator PLIP_INTRA_MOLECULAR_ANNOTATOR =
            new PLIPIntraMolecularAnnotator();

    public static void main(String[] args) {
        Map<Pair<AminoAcid.Family, AminoAcid.Family>, Double> map = new HashMap<>();
        Map<Pair<AminoAcid.Family, AminoAcid.Family>, Double> hydrogenbond = new HashMap<>();
        Map<Pair<AminoAcid.Family, AminoAcid.Family>, Double> hydrophobic = new HashMap<>();

        MembraneConstants.lines(DIRECTORY.resolve("ids.list"))
                .forEach(line -> {
                    System.out.println(line);
                    String pdbId = line.split("_")[0];
                    String chainId = line.split("_")[1];

                    Structure structure = StructureParser.source(DIRECTORY.resolve("pdb").resolve(pdbId + ".pdb"))
                            .minimalParsing(true)
                            .parse();
                    Chain chain = structure.select()
                            .chainName(chainId)
                            .asChain();

                    ORIENTATIONS_OF_PROTEINS_IN_MEMBRANES_ANNOTATOR.process(structure,
                            Jsoup.parse(MembraneConstants.lines(DIRECTORY.resolve("opm").resolve(pdbId + ".opm"))
                                    .collect(Collectors.joining(System.lineSeparator()))));
                    PLIP_INTRA_MOLECULAR_ANNOTATOR.process(chain,
                            Jsoup.parse(MembraneConstants.lines(DIRECTORY.resolve("plip").resolve(line + ".plip"))
                                    .collect(Collectors.joining(System.lineSeparator()))));

                    chain.aminoAcids()
                            .forEach(aminoAcid -> {
                                boolean tm = aminoAcid.getFeature(Topology.class).isTransmembrane();
                                if(!tm) {
                                    return;
                                }

                                PLIPInteractionContainer interactionContainer = aminoAcid.getFeature(PLIPInteractionContainer.class);
                                interactionContainer.getInteractions()
                                        .stream()
                                        .filter(MembraneConstants::isTmHelixInteraction)
                                        .forEach(interaction -> {
                                            AminoAcid.Family aa1 = AminoAcid.Family.resolveGroupPrototype(interaction.getPartner1().getGroupPrototype());
                                            AminoAcid.Family aa2 = AminoAcid.Family.resolveGroupPrototype(interaction.getPartner2().getGroupPrototype());

                                            Pair<AminoAcid.Family, AminoAcid.Family> pair1 = new Pair<>(aa1, aa2);
                                            Pair<AminoAcid.Family, AminoAcid.Family> pair2 = pair1.flip();

                                            // 0.5 to compensate for symmetry
                                            map.merge(pair1, 0.5, Double::sum);
                                            map.merge(pair2, 0.5, Double::sum);

                                            if(interaction instanceof HydrogenBond) {
                                                hydrogenbond.merge(pair1, 0.5, Double::sum);
                                                hydrogenbond.merge(pair2, 0.5, Double::sum);
                                            } else if(interaction instanceof HydrophobicInteraction) {
                                                hydrophobic.merge(pair1, 0.5, Double::sum);
                                                hydrophobic.merge(pair2, 0.5, Double::sum);
                                            }
                                        });
                            });
                });

        MembraneConstants.write(MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY.resolve("contact")
                        .resolve("aminoacids-contacts.mat"),
                composeOutput(map));
        MembraneConstants.write(MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY.resolve("contact")
                        .resolve("aminoacids-hydrogenbonds.mat"),
                composeOutput(hydrogenbond));
        MembraneConstants.write(MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY.resolve("contact")
                        .resolve("aminoacids-hydrophobic.mat"),
                composeOutput(hydrophobic));
    }

    private static String composeOutput(Map<Pair<AminoAcid.Family, AminoAcid.Family>, Double> map) {
        return "aa1;aa2;value" + System.lineSeparator() +
                AminoAcid.Family.canonicalAminoAcids()
                .flatMap(family1 -> AminoAcid.Family.canonicalAminoAcids()
                        .map(family2 -> new Pair<>(family1, family2)))
                .map(pair -> pair.getLeft().getOneLetterCode() + ";" +
                        pair.getRight().getOneLetterCode() + ";" +
                        StandardFormat.format(map.getOrDefault(pair, 0.0)))
                .collect(Collectors.joining(System.lineSeparator()));
    }
}
