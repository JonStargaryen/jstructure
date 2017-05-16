package design.aggregation;

import de.bioforscher.jstructure.feature.motif.SequenceMotifAnnotator;
import de.bioforscher.jstructure.feature.sse.SecStrucState;
import de.bioforscher.jstructure.feature.sse.SecondaryStructureAnnotator;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import de.bioforscher.jstructure.model.structure.family.GroupInformation;
import design.ProteinSource;
import design.parser.opm.OPMParser;
import design.parser.opm.TMHelix;

import java.io.IOException;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;

/**
 * Adds secondary structure information to proteins and prints all aggregated information to the console.
 * Moved to {@link ProteinSource}.
 * Created by S on 31.10.2016.
 */
@Deprecated
public class S04_ExtractSecondaryStructure {
    public static void main(String[] args) throws IOException {
        ProteinSource.loadProteins().forEach(S04_ExtractSecondaryStructure::printInformation);
    }

    public static void printInformation(Protein protein) {
        // print sequence
        System.out.println(Selection.on(protein)
                .aminoAcids()
                .asFilteredGroups()
                .map(Group::getGroupInformation)
                .map(GroupInformation::getOneLetterCode)
                .collect(Collectors.joining("", "seq: ", ""))
        );

        // print sequence motifs
        System.out.println(Selection.on(protein)
                .aminoAcids()
                .asFilteredGroups()
                .map(residue -> residue.getFeature(List.class,
                      SequenceMotifAnnotator.SEQUENCE_MOTIF))
                .map(entry -> Objects.isNull(entry) || entry.isEmpty() ? " " : "S")
                .collect(Collectors.joining("", "mot: " , ""))
        );

        // print topology
        System.out.println(Selection.on(protein)
                .aminoAcids()
                .asFilteredGroups()
                .map(residue -> residue.getFeature(TMHelix.class, OPMParser.TM_HELIX))
                .map(entry -> Objects.isNull(entry) ? " " : "M")
                .collect(Collectors.joining("", "top: ", ""))
        );

        // print secondary structure information
        System.out.println(Selection.on(protein)
                .aminoAcids()
                .asFilteredGroups()
                .map(residue -> residue.getFeature(SecStrucState.class,
                      SecondaryStructureAnnotator.SECONDARY_STRUCTURE_STATES))
                .map(state -> state.getSecondaryStructure().getReducedRepresentation())
                .collect(Collectors.joining("", "sse: ", ""))
        );
    }
}