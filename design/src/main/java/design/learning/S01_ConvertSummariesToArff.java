package design.learning;

import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import design.DesignConstants;

import java.nio.file.Path;
import java.nio.file.Paths;

/**
 * Are sequences connected to structure of motifs? GMLVQ the answer!
 * Created by bittrich on 12/19/16.
 */
public class S01_ConvertSummariesToArff {
    public static void main(String[] args) {
        for (String topology : DesignConstants.TOPOLOGIES) {
            for(SequenceMotifDefinition sequenceMotifDefinition : SequenceMotifDefinition.values()) {
                Path outPath = Paths.get(DesignConstants.ARFF_DIR + "sequence-to-clusters/" + topology + "/");
                DesignConstants.makeDirectoryIfAbsent(outPath);
                ClusterSummaryToArffConverter.writeToArff(Paths.get(DesignConstants.FRAGMENT_CLUSTERS + topology +
                        "/" + sequenceMotifDefinition.name() + "/" + DesignConstants.CLUSTER_SUMMARY),
                        Paths.get(outPath + "/" + sequenceMotifDefinition.name() + DesignConstants.ARFF_SUFFIX));
            }
        }
    }
}
