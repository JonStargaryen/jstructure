package design.aggregation;

import de.bioforscher.jstructure.feature.motif.SequenceMotif;
import de.bioforscher.jstructure.feature.motif.SequenceMotifAnnotator;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.Residue;
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
        // get proteins with sequence motif information
        loadProteins(false, true, true)
            .stream()
            .flatMap(Protein::residues)
            .filter(residue -> nonNull(residue.getFeature(List.class,
                    SequenceMotifAnnotator.FeatureNames.SEQUENCE_MOTIF)))
            .flatMap(residue -> residue.getFeature(List.class,
                    SequenceMotifAnnotator.FeatureNames.SEQUENCE_MOTIF).stream())
            .map(SequenceMotif.class::cast)
            .distinct()
            .forEach(S05_GroupByTopology::writeFragment);
    }

    private static void writeFragment(Object motifObject) {
        // why?
        SequenceMotif motif = (SequenceMotif) motifObject;
        System.out.println("writing file for " + motif);
        String proteinName = motif.getStartResidue().getParentChain().getParentProtein().getName();
        String topology = mapTopology(ProteinSource.determineTopologyGroup(motif));
        String chainId = motif.getStartResidue().getParentChain().getChainId();
        String content = motif.residues().map(Residue::composePDBRecord).collect(Collectors.joining(System.lineSeparator()));
        try {
            String filename = DesignConstants.MOTIF_FRAGMENT_BY_TOPOLOGY_DIR + topology + "/" + motif.getMotifDefinition() +
                    "-" + proteinName + "-" + chainId + "-" + motif.getStartResidue().getPdbName() +
                    motif.getStartResidue().getResidueNumber() + "-" + motif.getEndResidue().getPdbName() +
                    motif.getEndResidue().getResidueNumber() + DesignConstants.PDB_SUFFIX;
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
