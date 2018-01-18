package de.bioforscher.jstructure.feature.graphs;

import de.bioforscher.jstructure.mathematics.LinearAlgebra;
import de.bioforscher.jstructure.mathematics.graph.Graph;
import de.bioforscher.jstructure.model.feature.FeatureContainer;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

public class TopologicPropertyCalculator extends FeatureProvider {
    @Override
    protected void processInternally(Structure structure) {
        structure.chainsWithAminoAcids()
                .forEach(this::processInternally);
    }

    private void processInternally(Chain chain) {
        //TODO impl graph generation
        Graph<AminoAcid> plipGraph = null;
        Graph<AminoAcid> conventionalGraph = null;

        LinearAlgebra.GraphLinearAlgebra<AminoAcid> plipGraphCalculate = plipGraph.calculate();
        LinearAlgebra.GraphLinearAlgebra<AminoAcid> conventionalGraphCalculate = conventionalGraph.calculate();

        FeatureContainer chainFeatureContainer = chain.getFeatureContainer();
        chainFeatureContainer.addFeature(new PLIPChainTopologicProperties(this,
                plipGraphCalculate.graphPathLength(),
                plipGraphCalculate.clusteringCoefficient()));
        chainFeatureContainer.addFeature(new ConventionalChainTopologicProperties(this,
                conventionalGraphCalculate.graphPathLength(),
                conventionalGraphCalculate.clusteringCoefficient()));

        chain.aminoAcids()
                .forEach(aminoAcid -> processInternally(aminoAcid, plipGraphCalculate, conventionalGraphCalculate));
    }

    private void processInternally(AminoAcid aminoAcid,
                                   LinearAlgebra.GraphLinearAlgebra<AminoAcid> plipGraphCalculate,
                                   LinearAlgebra.GraphLinearAlgebra<AminoAcid> conventionalGraphCalculate) {
        FeatureContainer aminoAcidFeatureContainer = aminoAcid.getFeatureContainer();
        aminoAcidFeatureContainer.addFeature(new PLIPResidueTopologicProperties(this,
                plipGraphCalculate.betweenness(aminoAcid),
                plipGraphCalculate.closeness(aminoAcid),
                plipGraphCalculate.clusteringCoefficient(aminoAcid)));
        aminoAcidFeatureContainer.addFeature(new ConventionalResidueTopologicProperties(this,
                conventionalGraphCalculate.betweenness(aminoAcid),
                conventionalGraphCalculate.closeness(aminoAcid),
                conventionalGraphCalculate.clusteringCoefficient(aminoAcid)));
    }
}
