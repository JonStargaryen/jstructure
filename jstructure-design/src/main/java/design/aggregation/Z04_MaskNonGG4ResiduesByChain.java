package design.aggregation;

import de.bioforscher.jstructure.feature.motif.SequenceMotif;
import de.bioforscher.jstructure.feature.motif.SequenceMotifAnnotator;
import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.family.GroupInformation;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import de.bioforscher.jstructure.parser.ProteinParser;
import design.DesignConstants;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

/**
 * How do GG4 motifs interact with each other?
 * Created by S on 05.01.2017.
 */
public class Z04_MaskNonGG4ResiduesByChain {
    public static void main(String[] args) throws IOException {
        Files.lines(Paths.get(DesignConstants.NR_ALPHA_IDS))
                // select for non-comment lines
                .filter(DesignConstants.isCommentLine.negate())
                .forEach(line -> {
                    System.out.println(line);
                    String[] splitLine = line.split("_");
                    Protein protein = ProteinParser.source(splitLine[0]).parse();
                    // annotate sequence motifs in proteins
                    new SequenceMotifAnnotator().process(protein);

                    // mask residues not part sequence motif
                    Selection.on(protein)
                            .aminoAcids()
                            .asFilteredGroups()
                            .filter(Z04_MaskNonGG4ResiduesByChain::isNotPartOfGG4)
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
                        Files.write(Paths.get(DesignConstants.BASE_DIR + "masked-pdb-chains-gg4/" + protein.getName() + "_"
                                + chain.getChainId()+ DesignConstants.PDB_SUFFIX), chain.composePDBRecord().getBytes());
                    } catch(IOException e) {
                        e.printStackTrace();
                    }
                });
    }

    @SuppressWarnings("unchecked")
    private static boolean isNotPartOfGG4(Group residue) {
        List<SequenceMotif> motifAnnotations = residue.getFeature(List.class,
                SequenceMotifAnnotator.SEQUENCE_MOTIF);

        if(motifAnnotations == null || motifAnnotations.isEmpty()) {
            return true;
        }

        return motifAnnotations.stream()
                .filter(motif -> motif.getMotifDefinition().equals(SequenceMotifDefinition.GG4))
                .count() == 0;

    }
}
