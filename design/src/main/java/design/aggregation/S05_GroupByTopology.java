package design.aggregation;

import de.bioforscher.jstructure.feature.motif.SequenceMotif;
import de.bioforscher.jstructure.feature.motif.SequenceMotifAnnotator;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import design.DesignConstants;
import design.ProteinSource;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

import static design.ProteinSource.loadProteins;
import static java.util.Objects.nonNull;

/**
 * Separates sequence motifs by their topology (TM, nTM, transition: defined as a motif in more than 1 state). Extracted
 * motifs with coordinates are written to separate files.
 * Created by S on 01.11.2016.
 */
@SuppressWarnings("unchecked")
public class S05_GroupByTopology {
    public static void main(String[] args) throws IOException {
        // findAny proteins with sequence motif information
        loadProteins(false, true, true)
            .stream()
            .flatMap(protein -> Selection.on(protein)
                    .aminoAcids()
                    .filteredGroups())
            .filter(residue -> nonNull(residue.getFeature(List.class,
                    SequenceMotifAnnotator.FeatureNames.SEQUENCE_MOTIF)))
            .flatMap(residue -> residue.getFeature(List.class,
                    SequenceMotifAnnotator.FeatureNames.SEQUENCE_MOTIF).stream())
            .distinct()
            .map(SequenceMotif.class::cast)
            .forEach(S05_GroupByTopology::writeFragment);
    }

    private static void writeFragment(Object motifObject) {
        // why?
        SequenceMotif motif = (SequenceMotif) motifObject;
        System.out.println("writing file for " + motif);
        String proteinName = motif.getStartGroup().getParentChain().getParentProtein().getName();
        String topology = mapTopology(ProteinSource.determineTopologyGroup(motif));
        String chainId = motif.getStartGroup().getParentChain().getChainId();
        String content = motif.getGroupContainer()
                .groups()
                .map(Group::composePDBRecord)
                .collect(Collectors.joining(System.lineSeparator()));
        try {
            String filename = DesignConstants.MOTIF_FRAGMENT_BY_TOPOLOGY_DIR + topology + "/" + motif.getMotifDefinition() +
                    "-" + proteinName + "-" + chainId + "-" + motif.getStartGroup().getPdbName() +
                    motif.getStartGroup().getResidueNumber() + "-" + motif.getEndGroup().getPdbName() +
                    motif.getEndGroup().getResidueNumber() + DesignConstants.PDB_SUFFIX;
            Files.write(Paths.get(filename), content.getBytes());
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    private static String mapTopology(String key) {
        if(key.equals("o")) {
            return "ntm";
        }
        if(key.equals("I")) {
            return "tm";
        }
        return "trans";
    }
}
