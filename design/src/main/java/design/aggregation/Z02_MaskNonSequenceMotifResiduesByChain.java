package design.aggregation;

import de.bioforscher.jstructure.feature.motif.SequenceMotif;
import de.bioforscher.jstructure.feature.motif.SequenceMotifAnnotator;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import de.bioforscher.jstructure.model.structure.family.GroupInformation;
import de.bioforscher.jstructure.parser.ProteinParser;
import design.DesignConstants;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

/**
 * Masks all residues which are not part of a sequence motif by 'X'. Done to investigate structures further by Florian
 * and the itemset miner. This time getChain-specific.
 * Created by S on 07.11.2016.
 */
@Deprecated
public class Z02_MaskNonSequenceMotifResiduesByChain {
    public static void main(String[] args) throws IOException {
        Files.lines(Paths.get(DesignConstants.NR_ALPHA_IDS))
                // select for non-comment lines
                .filter(DesignConstants.isCommentLine.negate())
                .forEach(line -> {
                    System.out.println(line);
                    String[] splitLine = line.split("_");
                    Protein protein = ProteinParser.parseProteinById(splitLine[0]);
                    // annotate sequence motifs in proteins
                    new SequenceMotifAnnotator().process(protein);

                    // mask residues not part sequence motif
                    Selection.on(protein)
                            .aminoAcids()
                            .asFilteredGroups()
                            .filter(Z02_MaskNonSequenceMotifResiduesByChain::isNotStartOrEndOfSequenceMotif)
                            .forEach(residue -> residue.setGroupInformation(GroupInformation.UNKNOWN_AMINO_ACID));

                    System.out.printf("masked sequence for '%s' is:\n%s\n", protein.getName(), protein.getAminoAcidSequence());

                    // selectionBuilder appropriate getChain
                    String chainName;
                    try {
                        chainName = splitLine[1];
                    } catch (ArrayIndexOutOfBoundsException e) {
                        chainName = " ";
                    }
                    Chain chain = Selection.on(protein)
                            .chainName(chainName)
                            .asChain();

                    // write structures
                    try {
                        Files.write(Paths.get(DesignConstants.BASE_DIR + "masked-pdb-chains/" + protein.getName() + "_"
                                + chain.getChainId()+ DesignConstants.PDB_SUFFIX), chain.composePDBRecord().getBytes());
                    } catch(IOException e) {
                        e.printStackTrace();
                    }
                });
    }

    @SuppressWarnings("unchecked")
    private static boolean isNotStartOrEndOfSequenceMotif(Group residue) {
        List<SequenceMotif> motifAnnotations = residue.getFeature(List.class,
                SequenceMotifAnnotator.SEQUENCE_MOTIF);

        return motifAnnotations == null || motifAnnotations.isEmpty() || motifAnnotations.stream().filter(motif ->
                motif.getStartGroup().equals(residue) || motif.getEndGroup().equals(residue)).count() == 0;

    }
}
