package design.fragment.library;

import de.bioforscher.jstructure.alignment.AlignmentResult;
import de.bioforscher.jstructure.alignment.consensus.StructureCluster;
import de.bioforscher.jstructure.alignment.svd.SVDSuperimposer;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import design.DesignConstants;
import design.ProteinSource;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.List;

import static design.DesignConstants.DELIMITER;

/**
 * Previously, we approached the problem from a sequence perspective. However, why rely on Gerstein motifs, when all the
 * necessary information can be obtained directly from the data set.
 * TODO long-term residues surrounding the fragments (sequentially or spatially should be considered)
 * Created by S on 29.12.2016.
 */
public class S01_BuildLibrary {
    public static void main(String[] args) throws IOException {
        String basePath = DesignConstants.NAIVE_FRAGMENT_CLUSTERS + "/";
        DesignConstants.makeDirectoryIfAbsent(Paths.get(basePath));

        List<Protein> proteins = ProteinSource.loadProteins(false, false, true);
        List<StructureCluster> clusters = new StructureFragmentizer().fragmentize(proteins);

        // output consensus
        StringBuilder output = new StringBuilder();
        output.append("clusterId")
                .append(DELIMITER)
                .append("totalCount")
                .append(DELIMITER)
                .append("clusterSize")
                .append(DELIMITER)
                .append("fragmentId")
                .append(DELIMITER)
                .append("fragmentSequence")
                .append(DELIMITER)
                .append("rmsd")
                .append(System.lineSeparator());

        for (int i = 0; i < clusters.size(); i++) {
            if(i == 0 && clusters.get(0).getOriginalEntries().size() == 0) {
                // no rare clusters were merged - thus, no freak cluster at index 0
                continue;
            }

            String clusterPath = basePath + i + "/";
            DesignConstants.makeDirectoryIfAbsent(Paths.get(clusterPath));

            StructureCluster structureCluster = clusters.get(i);
            AtomContainer consensus = structureCluster.getConsensusRepresentation();
            DesignConstants.write(Paths.get(clusterPath + DesignConstants.CLUSTER_CONSENSUS),
                    consensus.composePDBRecord().getBytes());

            for(AtomContainer fragment : structureCluster.getOriginalEntries()) {
                // align fragment relative to consensus
                SVDSuperimposer svdSuperimposer = new SVDSuperimposer();
                AlignmentResult alignmentResult = svdSuperimposer.align(consensus, fragment);
                String rmsd = DesignConstants.DECIMAL_FORMAT.format(alignmentResult.getRmsd());

                // add to summary file
                output.append(i)
                        .append(DELIMITER)
                        .append(proteins.size())
                        .append(DELIMITER)
                        .append(structureCluster.getOriginalEntries().size())
                        .append(DELIMITER)
                        .append(fragment.getIdentifier())
                        .append(DELIMITER)
                        .append(((GroupContainer) fragment).getAminoAcidSequence())
                        .append(DELIMITER)
                        .append(rmsd)
                        .append(System.lineSeparator());

                // write aligned file
                DesignConstants.write(Paths.get(clusterPath + fragment.getIdentifier() + DesignConstants.PDB_SUFFIX),
                        alignmentResult.getAlignedQuery().composePDBRecord().getBytes());
            }
        }

        DesignConstants.write(Paths.get(basePath + DesignConstants.CLUSTER_SUMMARY), output.toString().getBytes());
    }
}
