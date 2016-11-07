package design.aggregation;

import de.bioforscher.jstructure.feature.motif.SequenceMotif;
import de.bioforscher.jstructure.feature.motif.SequenceMotifAnnotator;
import de.bioforscher.jstructure.model.structure.AminoAcid;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.Residue;
import de.bioforscher.jstructure.parser.ProteinParser;
import design.DesignConstants;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

/**
 * Masks all residues which are not part of a sequence motif by 'X'. Done to investigate structures further by Florian
 * and the itemset miner. This time chain-specific.
 * Created by S on 07.11.2016.
 */
public class Z02_MaskNonSequenceMotifResiduesByChain {
    public static void main(String[] args) throws IOException {
        Files.lines(Paths.get(DesignConstants.NR_ALPHA_IDS))
                // filter for non-comment lines
                .filter(DesignConstants.isCommentLine.negate())
                .forEach(line -> {
                    System.out.println(line);
                    String[] splitLine = line.split("_");
                    Protein protein = ProteinParser.parseProteinById(splitLine[0]);
                    // annotate sequence motifs in proteins
                    new SequenceMotifAnnotator().process(protein);

                    // mask residues not part sequence motif
                    protein.residues()
                            .filter(Z02_MaskNonSequenceMotifResiduesByChain::isNotStartOrEndOfSequenceMotif)
                            .forEach(residue -> residue.setAminoAcid(AminoAcid.UNKNOWN));

                    System.out.printf("masked sequence for '%s' is:\n%s\n", protein.getName(), protein.getSequence());

                    // select appropriate chain
                    Chain chain;
                    try {
                        chain = protein.chain(splitLine[1]).get();
                    } catch (ArrayIndexOutOfBoundsException e) {
                        chain = protein.chain(" ").get();
                    }

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
