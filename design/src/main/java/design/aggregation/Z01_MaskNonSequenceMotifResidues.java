package design.aggregation;

import de.bioforscher.jstructure.feature.motif.SequenceMotif;
import de.bioforscher.jstructure.feature.motif.SequenceMotifAnnotator;
import de.bioforscher.jstructure.model.structure.Group;
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
 * and the itemset miner.
 * Created by S on 29.10.2016.
 */
@Deprecated
public class Z01_MaskNonSequenceMotifResidues {
    public static void main(String[] args) throws IOException {
        Files.list(Paths.get(DesignConstants.PDB_DIR))
             // parse file
             .map(ProteinParser::parsePDBFile)
             //
             .forEach(protein -> {
                 // annotate sequence motifs in proteins
                 new SequenceMotifAnnotator().process(protein);

                 // mask residues not part sequence motif
                 Selection.on(protein)
                         .aminoAcids()
                         .asFilteredGroups()
                         .filter(Z01_MaskNonSequenceMotifResidues::isNotStartOrEndOfSequenceMotif)
                         .forEach(residue -> residue.setGroupInformation(GroupInformation.UNKNOWN_AMINO_ACID));

                 System.out.printf("masked sequence for '%s' is:\n%s\n", protein.getName(), protein.getAminoAcidSequence());

                 // write structures
                 try {
                     Files.write(Paths.get(DesignConstants.BASE_DIR + "masked-pdb/" + protein.getName() +
                             DesignConstants.PDB_SUFFIX), protein.composePDBRecord().getBytes());
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
