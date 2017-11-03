package de.bioforscher.jstructure.membrane.contact.definition;

import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.interactions.PLIPIntraMolecularAnnotator;
import de.bioforscher.jstructure.feature.topology.MembraneContainer;
import de.bioforscher.jstructure.feature.topology.OrientationsOfProteinsInMembranesAnnotator;
import de.bioforscher.jstructure.membrane.MembraneConstants;
import de.bioforscher.jstructure.membrane.GraphFactory;
import de.bioforscher.jstructure.mathematics.SetOperations;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.selection.IntegerRange;
import org.jsoup.Jsoup;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Visualize the impact of the contact definition on the way helices interact in 1fft.
 */
public class CreatePyMolFilesOfContactDefinitions {
    public static void main(String[] args) {
        String pdbId = "1fft";
        String chainId = "A";

        Structure structure = StructureParser.source(pdbId)
                .minimalParsing(true)
                .parse();
        Chain chain = structure.select()
                .chainName(chainId)
                .asChain();

        String opmDocument = MembraneConstants.lines(MembraneConstants.PDBTM_NR_ALPHA_DATASET_OPM_DIRECTORY.resolve(pdbId + ".opm"))
                .collect(Collectors.joining(System.lineSeparator()));
        new OrientationsOfProteinsInMembranesAnnotator().process(structure, Jsoup.parse(opmDocument));

        String plipDocument = MembraneConstants.lines(MembraneConstants.PDBTM_NR_ALPHA_DATASET_PLIP_DIRECTORY.resolve(pdbId + "_" + chainId + ".plip"))
                .collect(Collectors.joining(System.lineSeparator()));
        new PLIPIntraMolecularAnnotator().process(chain, Jsoup.parse(plipDocument));
        PLIPInteractionContainer plipInteractionContainer = chain.getFeature(PLIPInteractionContainer.class);

        List<IntegerRange> segments = structure.getFeature(MembraneContainer.class)
                .getTransMembraneSubunits()
                .stream()
                .filter(tms -> tms.getChainId().equals(chainId))
                .findFirst()
                .get()
                .getSegments();

        // determines the most interacting parts
        SetOperations.uniquePairsOf(segments)
                .forEach(segment -> {
                    List<Group> aminoAcids1 = chain.select()
                            .aminoAcids()
                            .residueNumber(segment.getLeft())
                            .asFilteredGroups()
                            .collect(Collectors.toList());
                    List<Group> aminoAcids2 = chain.select()
                            .aminoAcids()
                            .residueNumber(segment.getRight())
                            .asFilteredGroups()
                            .collect(Collectors.toList());

                    System.out.println(segment + " : " + SetOperations.cartesianProductOf(aminoAcids1, aminoAcids2)
                            .filter(pair -> plipInteractionContainer.areInContact(pair.getLeft(), pair.getRight()))
                            .count());
                });

        // (Range [414, 436], Range [451, 473]) : 7
        IntegerRange range1 = new IntegerRange(414, 436);
        IntegerRange range2 = new IntegerRange(451, 473);
        List<AminoAcid> aminoAcids1 = chain.select()
                .aminoAcids()
                .residueNumber(range1)
                .asFilteredGroups()
                .map(AminoAcid.class::cast)
                .collect(Collectors.toList());
        List<AminoAcid> aminoAcids2 = chain.select()
                .aminoAcids()
                .residueNumber(range2)
                .asFilteredGroups()
                .map(AminoAcid.class::cast)
                .collect(Collectors.toList());

        Stream.of(GraphFactory.InteractionScheme.values())
                .forEach(interactionScheme -> writePyMolFile(pdbId, chainId, range1, range2, aminoAcids1, aminoAcids2, interactionScheme));
    }

    private static void writePyMolFile(String pdbId,
                                       String chainId,
                                       IntegerRange range1,
                                       IntegerRange range2,
                                       List<AminoAcid> aminoAcids1,
                                       List<AminoAcid> aminoAcids2,
                                       GraphFactory.InteractionScheme interactionScheme) {
        String output = "fetch " + pdbId + ", async=0" + System.lineSeparator() +
                "bg_color white" + System.lineSeparator() +
                // hide non-relevant stuff
                "hide everything" + System.lineSeparator() +
                "select chain, chain " + chainId + System.lineSeparator() +
                "select residues1, chain " + chainId + " and resi " + range1.getLeft() + "-" + range1.getRight() + System.lineSeparator() +
                "select residues2, chain " + chainId + " and resi " + range2.getLeft() + "-" + range2.getRight() + System.lineSeparator() +
                "show_as cartoon, residues1" + System.lineSeparator() +
                "show_as cartoon, residues2" + System.lineSeparator() +
                SetOperations.cartesianProductOf(aminoAcids1, aminoAcids2)
                        .filter(pair -> interactionScheme.areInContact(pair.getLeft(), pair.getRight()))
                        .map(pair -> "distance chain " + chainId + " and resi " + pair.getLeft().getResidueIdentifier().getResidueNumber() + " and name CA, " +
                                "chain " + chainId + " and resi " + pair.getRight().getResidueIdentifier().getResidueNumber() + " and name CA")
                        .collect(Collectors.joining(System.lineSeparator())) + System.lineSeparator() +
                "set_view (\\" + System.lineSeparator() +
                     "0.660879016,   -0.730490386,    0.172111198,\\" + System.lineSeparator() +
                    "-0.749360025,   -0.629717827,    0.204723209,\\" + System.lineSeparator() +
                    "-0.041166950,   -0.264269710,   -0.963568509,\\" + System.lineSeparator() +
                     "0.000888713,    0.000804799,  -85.159309387,\\" + System.lineSeparator() +
                    "30.064170837,  291.794433594,  205.639144897,\\" + System.lineSeparator() +
                    "-1428.911987305, 1599.334594727,  -20.000000000)" + System.lineSeparator() +
                "color grey60" + System.lineSeparator() +
                "set ray_trace_mode, 3" + System.lineSeparator() +
                "hide labels";

        MembraneConstants.write(MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY.resolve("pymol").resolve(pdbId + "_"
                + chainId + "_" + interactionScheme + ".pml"), output);
    }
}
