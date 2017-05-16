package aars;

import de.bioforscher.jstructure.parser.pdb.PDBDatabaseQuery;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Ensure the single-linkage clustering is valid. To achieve that compare our clusters to that of the PDB.
 * Created by bittrich on 5/8/17.
 */
public class CheckIntegrityOfClustering {
    private static List<String> allIds;
    private static List<List<String>> clusters;
    private static Map<String, List<String>> pdbClusters;

    public static void main(String[] args) throws IOException {
        // load clusters of single-linkage clustering
        clusters = Files.lines(Paths.get(System.getProperty("user.home") + "/git/aars_analysis/data/sequence_clusters_cleansed.dat"))
                .filter(line -> !line.startsWith("#"))
                .map(line -> line.replace("[", "").replace("]", ""))
                .map(line -> line.split(", "))
                .map(split -> Stream.of(split).collect(Collectors.toList()))
                .collect(Collectors.toList());

        // extract all unique chainIds
        allIds = clusters.stream()
                .flatMap(Collection::stream)
                .collect(Collectors.toList());

        // fetch all pdb clusters
        pdbClusters = allIds.stream()
                .collect(Collectors.toMap(Function.identity(), id -> PDBDatabaseQuery.fetchSequenceCluster(id.split("_")[0], id.split("_")[1], 95)));

        // retain but considered structures, drop potential additional ballast added by the pdb query
        pdbClusters.values().forEach(list -> list.retainAll(allIds));

        clusters.forEach(CheckIntegrityOfClustering::checkIntegrity);
    }

    private static void checkIntegrity(List<String> cluster) {
        if(dissimilar(cluster, pdbClusters.get(cluster.get(0)))) {
            System.out.println("slc: " + cluster);
            for (String id : cluster) {
                System.out.println("pdb: " + pdbClusters.get(id));
            }
        }
    }

    private static boolean dissimilar(List<String> slCluster, List<String> pdbCluster) {
        if(slCluster.size() != pdbCluster.size()) {
            return true;
        }

        for(String id : slCluster) {
            if(!pdbCluster.contains(id)) {
                return true;
            }
        }

        for(String id : pdbCluster) {
            if(!slCluster.contains(id)) {
                return true;
            }
        }

        return false;
    }
}
