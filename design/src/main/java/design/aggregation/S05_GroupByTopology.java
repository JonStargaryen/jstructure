package design.aggregation;

import de.bioforscher.jstructure.feature.motif.SequenceMotif;
import de.bioforscher.jstructure.feature.motif.SequenceMotifAnnotator;
import de.bioforscher.jstructure.model.structure.Group;
import design.DesignConstants;
import design.ProteinSource;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.Collection;
import java.util.stream.Collectors;

import static design.ProteinSource.loadProteins;

/**
 * Separates sequence motifs by their topology (TM, nTM, transition: defined as a motif in more than 1 state). Extracted
 * motifs with coordinates are written to separate files.
 *
 * Diary:
 * <pre>
 *     Writes extracted sequence motif fragments to a dedicated directory.
 *
 *     12.14.16
 *      - rerun
 *      - added check to validate fragments (however, the check never fails as only the {@link SequenceMotifAnnotator} was bugged)
 *      - observations are now centered
 * </pre>
 * Created by S on 01.11.2016.
 */
@SuppressWarnings("unchecked")
public class S05_GroupByTopology {
    public static void main(String[] args) throws IOException {
        // findAny proteins with sequence motif information
        loadProteins(false, true, true)
                .stream()
                .map(protein -> protein.getFeatureAsList(SequenceMotif.class, SequenceMotifAnnotator.FeatureNames.SEQUENCE_MOTIF))
                .flatMap(Collection::stream)
                .forEach(S05_GroupByTopology::writeFragment);
    }

    private static void writeFragment(Object motifObject) {
        // why?
        SequenceMotif motif = (SequenceMotif) motifObject;
        System.out.println("writing file for " + motif);
        String proteinName = motif.getStartGroup().getParentChain().getParentProtein().getName();
        String topology = mapTopology(ProteinSource.determineTopologyGroup(motif));
        String chainId = motif.getStartGroup().getParentChain().getChainId();

        // center atoms
//        CoordinateManipulations.center(motif.getGroupContainer());

        String content = motif.getGroupContainer()
                .groups()
                .map(Group::composePDBRecord)
                .collect(Collectors.joining(""));

        // this check is not needed anymore
//        if(isErroneous(motif)) {
//            return;
//        }

        String filename = DesignConstants.MOTIF_FRAGMENT_BY_TOPOLOGY_DIR + topology + "/" + motif.getMotifDefinition() +
                "-" + proteinName + "-" + chainId + "-" + motif.getStartGroup().getPdbName() +
                motif.getStartGroup().getResidueNumber() + "-" + motif.getEndGroup().getPdbName() +
                motif.getEndGroup().getResidueNumber() + DesignConstants.PDB_SUFFIX;
        DesignConstants.write(Paths.get(filename), content.getBytes());
    }

    private static boolean isErroneous(SequenceMotif sequenceMotif) {
        Group startResidue = sequenceMotif.getStartGroup();
        Group endResidue = sequenceMotif.getEndGroup();
        int expectedLength = Integer.valueOf(sequenceMotif.getMotifDefinition().name().substring(2));

        if(!startResidue.getParentChain().equals(endResidue.getParentChain())) {
            System.out.println("parent chains are not identical");
            return true;
        }

        if(endResidue.getResidueNumber() - startResidue.getResidueNumber() != expectedLength) {
            System.out.println("residue numbers do not match");
            return true;
        }

        if(sequenceMotif.getGroupContainer().getGroups().size() != expectedLength + 1) {
            System.out.println("number of groups too small");
            return true;
        }

        if(sequenceMotif.getGroupContainer().getAtoms().size() == 0) {
            System.out.println("no ATOM content");
            return true;
        }

        return false;
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
