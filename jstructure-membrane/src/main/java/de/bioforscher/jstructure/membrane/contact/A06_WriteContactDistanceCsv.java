package de.bioforscher.jstructure.membrane.contact;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.feature.interactions.PLIPInteraction;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.interactions.PLIPIntraMolecularAnnotator;
import de.bioforscher.jstructure.feature.rigidity.BackboneRigidity;
import de.bioforscher.jstructure.feature.rigidity.DynaMineBridge;
import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.feature.topology.MembraneContainer;
import de.bioforscher.jstructure.feature.topology.OrientationsOfProteinsInMembranesAnnotator;
import de.bioforscher.jstructure.membrane.MembraneConstants;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.selection.IntegerRange;
import org.jsoup.Jsoup;

import java.nio.file.Path;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

/**
 * Do interactions types differ in their distance? Well, they should for sure.
 */
public class A06_WriteContactDistanceCsv {
    private static final OrientationsOfProteinsInMembranesAnnotator ORIENTATIONS_OF_PROTEINS_IN_MEMBRANES_ANNOTATOR =
            new OrientationsOfProteinsInMembranesAnnotator();
    private static final PLIPIntraMolecularAnnotator PLIP_INTRA_MOLECULAR_ANNOTATOR =
            new PLIPIntraMolecularAnnotator();
    private static final DynaMineBridge DYNA_MINE_BRIDGE =
            new DynaMineBridge();

    public static void main(String[] args) {
        String header = "pdbId;chainId;res1;aa1;res2;aa2;topology;size1;size2;interaction;distance;rigidity1;rigidity2" + System.lineSeparator();

        String alpha = composeCsv(MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY);
        MembraneConstants.write(MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY.resolve("contact").resolve("contact-distance.csv"),
                header + alpha);

        String beta = composeCsv(MembraneConstants.PDBTM_NR_BETA_DATASET_DIRECTORY);
        MembraneConstants.write(MembraneConstants.PDBTM_NR_BETA_DATASET_DIRECTORY.resolve("contact").resolve("contact-distance.csv"),
                header + beta);
    }

    private static String composeCsv(Path directory) {
        return MembraneConstants.lines(directory.resolve("ids.list"))
                .map(line -> handleLine(line, directory))
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.joining(System.lineSeparator()));
    }

    private static Optional<String> handleLine(String line, Path directory) {
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
            DYNA_MINE_BRIDGE.process(chain,
                    MembraneConstants.lines(directory.resolve("dynamine").resolve(line + "_backbone.pred"))
                            .collect(Collectors.joining(System.lineSeparator())));

            MembraneContainer membraneContainer = structure.getFeature(MembraneContainer.class);
            List<PLIPInteraction> interactions = chain.getFeature(PLIPInteractionContainer.class).getInteractions();

            double min = chain.aminoAcids()
                    .mapToDouble(aa -> aa.getFeature(BackboneRigidity.class).getBackboneRigidity())
                    .min()
                    .getAsDouble();
            double max = chain.aminoAcids()
                    .mapToDouble(aa -> aa.getFeature(BackboneRigidity.class).getBackboneRigidity())
                    .max()
                    .getAsDouble();

            return Optional.of(interactions.stream()
                    .map(interaction -> handleInteraction(interaction, membraneContainer, pdbId, chainId, min, max))
                    .filter(Optional::isPresent)
                    .map(Optional::get)
                    .collect(Collectors.joining(System.lineSeparator())));
        } catch (Exception e) {
            e.printStackTrace();
            return Optional.empty();
        }
    }

    private static Optional<String> handleInteraction(PLIPInteraction interaction,
                                                      MembraneContainer membraneContainer,
                                                      String pdbId,
                                                      String chainId,
                                                      double min,
                                                      double max) {
        Group partner1 = interaction.getPartner1();
        Group partner2 = interaction.getPartner2();
        Optional<IntegerRange> segment1 = membraneContainer.getEmbeddingTransmembraneSegment(partner1);
        Optional<IntegerRange> segment2 = membraneContainer.getEmbeddingTransmembraneSegment(partner2);
        String type;
        int size1;
        int size2;

        // ensure both are amino acids
        if(!(partner1 instanceof AminoAcid) || !(partner2 instanceof AminoAcid)) {
            return Optional.empty();
        }

        AminoAcid aminoAcid1 = (AminoAcid) partner1;
        AminoAcid aminoAcid2 = (AminoAcid) partner2;

        GenericSecondaryStructure.SecondaryStructureElement secondaryStructure1 = partner1.getFeature(GenericSecondaryStructure.class)
                .getSurroundingSecondaryStructureElement(aminoAcid1);
        GenericSecondaryStructure.SecondaryStructureElement secondaryStructure2 = partner2.getFeature(GenericSecondaryStructure.class)
                .getSurroundingSecondaryStructureElement(aminoAcid2);

        // ensure both are transmembrane helices and not in the same segment
        if(segment1.isPresent() &&  segment2.isPresent() && !segment1.equals(segment2)) {
            type = "I";
            size1 = segment1.get().getSize();
            size2 = segment2.get().getSize();
        } else if(secondaryStructure1.getReducedType().equals("H") && secondaryStructure2.getReducedType().equals("H") && !secondaryStructure1.equals(secondaryStructure2)) {
            type ="o";
            size1 = secondaryStructure1.getSize();
            size2 = secondaryStructure2.getSize();
        } else {
            return Optional.empty();
        }

        double caDistance = aminoAcid1.getCa()
                .calculate()
                .distance(aminoAcid2.getCa());
        double rigidity1 = minMaxNormalize(aminoAcid1.getFeature(BackboneRigidity.class).getBackboneRigidity(), min, max);
        double rigidity2 = minMaxNormalize(aminoAcid2.getFeature(BackboneRigidity.class).getBackboneRigidity(), min, max);

        String output = pdbId + ";" +
                chainId + ";" +
                partner1.getResidueIdentifier() + ";" +
                aminoAcid1.getOneLetterCode() + ";" +
                partner2.getResidueIdentifier() + ";" +
                aminoAcid2.getOneLetterCode() + ";" +
                type + ";" +
                size1 + ";" +
                size2 + ";" +
                interaction.getClass().getSimpleName() + ";" +
                StandardFormat.format(caDistance) + ";" +
                StandardFormat.format(rigidity1) + ";" +
                StandardFormat.format(rigidity2);

        return Optional.of(output);
    }

    private static double minMaxNormalize(double v, double min, double max) {
        return (v - min) / (max  - min);
    }
}
