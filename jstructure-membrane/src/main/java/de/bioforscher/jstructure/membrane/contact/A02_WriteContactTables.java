package de.bioforscher.jstructure.membrane.contact;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.feature.interactions.HydrogenBond;
import de.bioforscher.jstructure.feature.interactions.HydrophobicInteraction;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.interactions.PLIPIntraMolecularAnnotator;
import de.bioforscher.jstructure.feature.topology.OrientationsOfProteinsInMembranesAnnotator;
import de.bioforscher.jstructure.feature.topology.Topology;
import de.bioforscher.jstructure.mathematics.Pair;
import de.bioforscher.jstructure.mathematics.SetOperations;
import de.bioforscher.jstructure.membrane.MembraneConstants;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.jsoup.Jsoup;

import java.nio.file.Path;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public class A02_WriteContactTables {
    private static final Path DIRECTORY = MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY;
    private static final OrientationsOfProteinsInMembranesAnnotator ORIENTATIONS_OF_PROTEINS_IN_MEMBRANES_ANNOTATOR =
            new OrientationsOfProteinsInMembranesAnnotator();
    private static final PLIPIntraMolecularAnnotator PLIP_INTRA_MOLECULAR_ANNOTATOR =
            new PLIPIntraMolecularAnnotator();

    public static void main(String[] args) {
        Map<Pair<AminoAcid.Family, AminoAcid.Family>, Integer> map = new HashMap<>();
        Map<Pair<AminoAcid.Family, AminoAcid.Family>, Integer> hydrogenbond = new HashMap<>();
        Map<Pair<AminoAcid.Family, AminoAcid.Family>, Integer> hydrophobic = new HashMap<>();
        Map<AminoAcid.Family, Integer> occurrence = new HashMap<>();

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

                                // normalized by frequency of both interaction partners in the TM layer
                                occurrence.merge(AminoAcid.Family.resolveGroupPrototype(aminoAcid.getGroupPrototype()), 1, Integer::sum);

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
                                            map.merge(pair1, 1, Integer::sum);
                                            map.merge(pair2, 1, Integer::sum);

                                            if(interaction instanceof HydrogenBond) {
                                                hydrogenbond.merge(pair1, 1, Integer::sum);
                                                hydrogenbond.merge(pair2, 1, Integer::sum);
                                            } else if(interaction instanceof HydrophobicInteraction) {
                                                hydrophobic.merge(pair1, 1, Integer::sum);
                                                hydrophobic.merge(pair2, 1, Integer::sum);
                                            }
                                        });
                            });
                });

        MembraneConstants.write(MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY.resolve("results")
                        .resolve("contacts-freq.tab"),
                composeFrequencyOutput(map, occurrence));
        MembraneConstants.write(MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY.resolve("results")
                        .resolve("hydrogenbonds-freq.tab"),
                composeFrequencyOutput(hydrogenbond, occurrence));
        MembraneConstants.write(MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY.resolve("results")
                        .resolve("hydrophobic-freq.tab"),
                composeFrequencyOutput(hydrophobic, occurrence));
    }

    private static String composeFrequencyOutput(Map<Pair<AminoAcid.Family, AminoAcid.Family>, Integer> map,
                                                 Map<AminoAcid.Family, Integer> occurrence) {
        List<AminoAcid.Family> aminoAcids = AminoAcid.Family.canonicalAminoAcids().collect(Collectors.toList());
        return "aa1;aa2;abs;norm" + System.lineSeparator() +
                SetOperations.cartesianProductOf(aminoAcids, aminoAcids)
                .filter(pair -> pair.getRight().getOneLetterCode().compareTo(pair.getLeft().getOneLetterCode()) >= 0)
                .map(pair -> pair.getLeft().getOneLetterCode() + ";" + pair.getRight().getOneLetterCode() + ";" +
                        StandardFormat.format(0.5 * 100 * 100 * 100 * (map.getOrDefault(pair, 0))) + ";" +
                        StandardFormat.format(0.5 * 100 * 100 * 100 * (map.getOrDefault(pair, 0) /
                                ((double) occurrence.getOrDefault(pair.getLeft(), 1) * (double) occurrence.getOrDefault(pair.getRight(), 1)))))
                .collect(Collectors.joining(System.lineSeparator()));
    }
}
