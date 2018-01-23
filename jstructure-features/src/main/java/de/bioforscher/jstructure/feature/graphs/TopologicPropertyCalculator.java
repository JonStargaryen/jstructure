package de.bioforscher.jstructure.feature.graphs;

import de.bioforscher.jstructure.mathematics.LinearAlgebra;
import de.bioforscher.jstructure.mathematics.graph.Graph;
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

//    private void processInternally(Chain chain) {
//        Graph<AminoAcid> fullPlipGraph = ProteinGraphFactory.createProteinGraph(chain,
//                ProteinGraphFactory.InteractionScheme.SALENTIN2015);
//        Graph<AminoAcid> hydrogenPlipGraph = ProteinGraphFactory.createProteinGraph(chain,
//                ProteinGraphFactory.InteractionScheme.SALENTIN2015_HYDROGEN_BONDS);
//        Graph<AminoAcid> hydrophobicPlipGraph = ProteinGraphFactory.createProteinGraph(chain,
//                ProteinGraphFactory.InteractionScheme.SALENTIN2015_HYDROPHOBIC_INTERACTION);
//        Graph<AminoAcid> conventionalGraph = ProteinGraphFactory.createProteinGraph(chain,
//                ProteinGraphFactory.InteractionScheme.CALPHA8);
//
//        LinearAlgebra.GraphLinearAlgebra<AminoAcid> fullPlipGraphCalculate = fullPlipGraph.calculate();
//        LinearAlgebra.GraphLinearAlgebra<AminoAcid> hydrogenPlipGraphCalculate = hydrogenPlipGraph.calculate();
//        LinearAlgebra.GraphLinearAlgebra<AminoAcid> hydrophobicPlipGraphCalculate = hydrophobicPlipGraph.calculate();
//        LinearAlgebra.GraphLinearAlgebra<AminoAcid> conventionalGraphCalculate = conventionalGraph.calculate();
//
//        ChainTopologicProperties fullPlipChainTopologicProperties = new ChainTopologicProperties(this,
//                fullPlipGraphCalculate.averageGraphPathLength(),
//                fullPlipGraphCalculate.clusteringCoefficient());
//        ChainTopologicProperties hydrogenPlipChainTopologicProperties = new ChainTopologicProperties(this,
//                hydrogenPlipGraphCalculate.averageGraphPathLength(),
//                hydrogenPlipGraphCalculate.clusteringCoefficient());
//        ChainTopologicProperties hydrophobicPlipChainTopologicProperties = new ChainTopologicProperties(this,
//                hydrophobicPlipGraphCalculate.averageGraphPathLength(),
//                hydrophobicPlipGraphCalculate.clusteringCoefficient());
//        ChainTopologicProperties conventionalChainTopologicProperties = new ChainTopologicProperties(this,
//                conventionalGraphCalculate.averageGraphPathLength(),
//                conventionalGraphCalculate.clusteringCoefficient());
//
//        chain.getFeatureContainer().addFeature(new ChainTopologicPropertiesContainer(this,
//                fullPlipChainTopologicProperties,
//                hydrogenPlipChainTopologicProperties,
//                hydrophobicPlipChainTopologicProperties,
//                conventionalChainTopologicProperties));
//
//        chain.aminoAcids()
//                .forEach(aminoAcid -> processInternally(aminoAcid,
//                        fullPlipGraphCalculate,
//                        hydrogenPlipGraphCalculate,
//                        hydrophobicPlipGraphCalculate,
//                        conventionalGraphCalculate));
//    }
//
//    private void processInternally(AminoAcid aminoAcid,
//                                   LinearAlgebra.GraphLinearAlgebra<AminoAcid> fullPlipGraphCalculate,
//                                   LinearAlgebra.GraphLinearAlgebra<AminoAcid> hydrogenPlipGraphCalculate,
//                                   LinearAlgebra.GraphLinearAlgebra<AminoAcid> hydrophobicPlipGraphCalculate,
//                                   LinearAlgebra.GraphLinearAlgebra<AminoAcid> conventionalGraphCalculate) {
//        ResidueTopologicProperties fullPlipResidueTopologicProperties = new ResidueTopologicProperties(this,
//                fullPlipGraphCalculate.betweenness(aminoAcid),
//                fullPlipGraphCalculate.closeness(aminoAcid),
//                fullPlipGraphCalculate.clusteringCoefficient(aminoAcid));
//        ResidueTopologicProperties hydrogenPlipResidueTopologicProperties = new ResidueTopologicProperties(this,
//                hydrogenPlipGraphCalculate.betweenness(aminoAcid),
//                hydrogenPlipGraphCalculate.closeness(aminoAcid),
//                hydrogenPlipGraphCalculate.clusteringCoefficient(aminoAcid));
//        ResidueTopologicProperties hydrophobicPlipResidueTopologicProperties = new ResidueTopologicProperties(this,
//                hydrophobicPlipGraphCalculate.betweenness(aminoAcid),
//                hydrophobicPlipGraphCalculate.closeness(aminoAcid),
//                hydrophobicPlipGraphCalculate.clusteringCoefficient(aminoAcid));
//        ResidueTopologicProperties conventionalResidueTopologicProperties = new ResidueTopologicProperties(this,
//                conventionalGraphCalculate.betweenness(aminoAcid),
//                conventionalGraphCalculate.closeness(aminoAcid),
//                conventionalGraphCalculate.clusteringCoefficient(aminoAcid));
//        aminoAcid.getFeatureContainer().addFeature(new ResidueTopologicPropertiesContainer(this,
//                fullPlipResidueTopologicProperties,
//                hydrogenPlipResidueTopologicProperties,
//                hydrophobicPlipResidueTopologicProperties,
//                conventionalResidueTopologicProperties));
//    }

    private void processInternally(Chain chain) {
        Graph<AminoAcid> conventionalGraph = ProteinGraphFactory.createProteinGraph(chain,
                ProteinGraphFactory.InteractionScheme.CALPHA8);

        LinearAlgebra.GraphLinearAlgebra<AminoAcid> conventionalGraphCalculate = conventionalGraph.calculate();

        ChainTopologicProperties conventionalChainTopologicProperties = new ChainTopologicProperties(this,
                conventionalGraphCalculate.averageGraphPathLength(),
                conventionalGraphCalculate.clusteringCoefficient());

        chain.getFeatureContainer().addFeature(new ChainTopologicPropertiesContainer(this,
                null,
                null,
                null,
                conventionalChainTopologicProperties));

        chain.aminoAcids()
                .forEach(aminoAcid -> processInternally(aminoAcid,
                        conventionalGraphCalculate));
    }

    private void processInternally(AminoAcid aminoAcid,
                                   LinearAlgebra.GraphLinearAlgebra<AminoAcid> conventionalGraphCalculate) {
        ResidueTopologicProperties conventionalResidueTopologicProperties = new ResidueTopologicProperties(this,
                conventionalGraphCalculate.betweenness(aminoAcid),
                conventionalGraphCalculate.closeness(aminoAcid),
                conventionalGraphCalculate.clusteringCoefficient(aminoAcid));
        aminoAcid.getFeatureContainer().addFeature(new ResidueTopologicPropertiesContainer(this,
                null,
                null,
                null,
                conventionalResidueTopologicProperties));
    }
}
