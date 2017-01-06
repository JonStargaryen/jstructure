package design.visualization;

import de.bioforscher.jstructure.feature.motif.SequenceMotif;
import de.bioforscher.jstructure.feature.motif.SequenceMotifAnnotator;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import de.bioforscher.jstructure.model.structure.family.GroupInformation;
import de.bioforscher.jstructure.parser.ProteinParser;
import design.DesignConstants;

import java.nio.file.Paths;
import java.util.List;

/**
 * Take an OPM structure with annotated membrane and visualize all sequence motifs.
 * Created by bittrich on 12/19/16.
 */
public class S04_MaskSequencePositions {
    public static void main(String[] args) {
        Protein protein = ProteinParser.parsePDBFile("/home/bittrich/Downloads/1m0l.pdb");

        // annotate sequence motifs in proteins
        new SequenceMotifAnnotator().process(protein);

        // mask residues not part sequence motif
        Selection.on(protein)
                .aminoAcids()
                .asFilteredGroups()
                .filter(S04_MaskSequencePositions::isNotStartOrEndOfSequenceMotif)
                .forEach(residue -> residue.setGroupInformation(GroupInformation.UNKNOWN_AMINO_ACID));

        System.out.printf("masked sequence for '%s' is:\n%s\n", protein.getName(), protein.getAminoAcidSequence());

        // write structures
        DesignConstants.write(Paths.get("/home/bittrich/Downloads/1m0l-masked.pdb"), protein.composePDBRecord().getBytes());
    }

    @SuppressWarnings("unchecked")
    private static boolean isNotStartOrEndOfSequenceMotif(Group residue) {
        List<SequenceMotif> motifAnnotations = residue.getFeature(List.class,
                SequenceMotifAnnotator.FeatureNames.SEQUENCE_MOTIF);

        return motifAnnotations == null || motifAnnotations.isEmpty() || motifAnnotations.stream().filter(motif ->
                motif.getStartGroup().equals(residue) || motif.getEndGroup().equals(residue)).count() == 0;

    }
}
