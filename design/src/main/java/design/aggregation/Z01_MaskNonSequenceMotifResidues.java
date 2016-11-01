package design.aggregation;

import de.bioforscher.jstructure.feature.motif.SequenceMotif;
import de.bioforscher.jstructure.feature.motif.SequenceMotifAnnotator;
import de.bioforscher.jstructure.model.structure.AminoAcid;
import de.bioforscher.jstructure.model.structure.Residue;
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
                 protein.residues()
                        .filter(Z01_MaskNonSequenceMotifResidues::isNotStartOrEndOfSequenceMotif)
                        .forEach(residue -> residue.setAminoAcid(AminoAcid.UNKNOWN));

                 System.out.printf("masked sequence for '%s' is:\n%s\n", protein.getName(), protein.getSequence());

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
    private static boolean isNotStartOrEndOfSequenceMotif(Residue residue) {
        List<SequenceMotif> motifAnnotations = residue.getFeature(List.class,
                SequenceMotifAnnotator.FeatureNames.SEQUENCE_MOTIF);

        if(motifAnnotations == null || motifAnnotations.isEmpty()) {
            return true;
        }

        return motifAnnotations.stream()
                               .filter(motif -> motif.getStartResidue().equals(residue) ||
                                                motif.getEndResidue().equals(residue))
                               .count() == 0;
    }
}
