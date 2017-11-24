package de.bioforscher.jstructure.contacts.algorithms;

import de.bioforscher.jstructure.mathematics.Pair;
import de.bioforscher.jstructure.mathematics.graph.Edge;
import de.bioforscher.jstructure.mathematics.graph.Graph;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.util.Comparator;
import java.util.List;
import java.util.stream.Stream;

/**
 * Identifies structural essence of a contact map by removing redundant edges according to a common neighborhood
 * criterion.
 *
 * Sathyapriya R, Duarte JM, Stehr H, Filippis I, Lappe M (2009) Defining an Essence of Structure Determining Residue
 * Contacts in Proteins. PLoS Comput Biol 5(12): e1000584. doi:10.1371/journal.pcbi.1000584
 */
public class ConePeeling {
    public Graph<AminoAcid> conePeelProteinGraph(Graph<AminoAcid> originalGraph) {
        Graph<AminoAcid> graph = new Graph<>(originalGraph);

        // sort nodes in descending order by degree
        List<AminoAcid> nodes = graph.getNodes();
        nodes.sort(Comparator.comparingInt(graph::getDegreeOf));

        // for all nodes: get edges incident on it
        for(AminoAcid node : nodes) {
            List<Edge<AminoAcid>> incidentEdges = graph.getEdgesFor(node);
            // sort by common neighborhood
            incidentEdges.sort(Comparator.comparingInt(graph::getSizeOfCommonNeighborhood));

            incidentEdges.stream()
                    .flatMap(graph::commonNeighborhood)
                    .filter(this::neighborhoodToRemove)
                    .flatMap(pair -> Stream.of(pair.getLeft(), pair.getRight()))
                    .forEach(graph::remove);
        }

        return graph;
    }

    private boolean neighborhoodToRemove(Pair<Edge<AminoAcid>, Edge<AminoAcid>> pair) {
        AminoAcid central = pair.getLeft().getLeft();
        AminoAcid aminoAcid1 = pair.getLeft().getRight();
        AminoAcid aminoAcid2 = pair.getRight().getRight();
        return getSequenceRange(central, aminoAcid1) <= 3 || getSequenceRange(central, aminoAcid2) <= 3;
    }

    private int getSequenceRange(AminoAcid aminoAcid1, AminoAcid aminoAcid2) {
        return Math.abs(aminoAcid1.getResidueIndex() - aminoAcid2.getResidueIndex());
    }
}
