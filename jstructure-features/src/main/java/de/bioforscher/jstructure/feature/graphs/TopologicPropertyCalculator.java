package de.bioforscher.jstructure.feature.graphs;

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
        ProteinGraph fullPlipGraph = ProteinGraphFactory.createProteinGraph(chain,
                ProteinGraphFactory.InteractionScheme.SALENTIN2015);
        ProteinGraph hydrogenPlipGraph = ProteinGraphFactory.createProteinGraph(chain,
                ProteinGraphFactory.InteractionScheme.SALENTIN2015_HYDROGEN_BONDS);
        ProteinGraph hydrophobicPlipGraph = ProteinGraphFactory.createProteinGraph(chain,
                ProteinGraphFactory.InteractionScheme.SALENTIN2015_HYDROPHOBIC_INTERACTION);
        ProteinGraph conventionalGraph = ProteinGraphFactory.createProteinGraph(chain,
                ProteinGraphFactory.InteractionScheme.CALPHA8);

        ProteinGraphCalculations fullPlipCalculations = new ProteinGraphCalculations(fullPlipGraph);
        ProteinGraphCalculations hydrogenPlipCalculations = new ProteinGraphCalculations(hydrogenPlipGraph);
        ProteinGraphCalculations hydrophobicPlipCalculations = new ProteinGraphCalculations(hydrophobicPlipGraph);
        ProteinGraphCalculations conventionalCalculations = new ProteinGraphCalculations(conventionalGraph);

        ChainTopologicProperties fullPlipChainTopologicProperties = new ChainTopologicProperties(this,
                fullPlipCalculations.averageGraphPathLength(),
                fullPlipCalculations.clusteringCoefficient());
        ChainTopologicProperties hydrogenPlipChainTopologicProperties = new ChainTopologicProperties(this,
                hydrogenPlipCalculations.averageGraphPathLength(),
                hydrogenPlipCalculations.clusteringCoefficient());
        ChainTopologicProperties hydrophobicPlipChainTopologicProperties = new ChainTopologicProperties(this,
                hydrophobicPlipCalculations.averageGraphPathLength(),
                hydrophobicPlipCalculations.clusteringCoefficient());
        ChainTopologicProperties conventionalChainTopologicProperties = new ChainTopologicProperties(this,
                conventionalCalculations.averageGraphPathLength(),
                conventionalCalculations.clusteringCoefficient());

        chain.getFeatureContainer().addFeature(new ChainTopologicPropertiesContainer(this,
                fullPlipChainTopologicProperties,
                hydrogenPlipChainTopologicProperties,
                hydrophobicPlipChainTopologicProperties,
                conventionalChainTopologicProperties));

        chain.aminoAcids()
                .forEach(aminoAcid -> processInternally(aminoAcid,
                        fullPlipCalculations,
                        hydrogenPlipCalculations,
                        hydrophobicPlipCalculations,
                        conventionalCalculations));
    }

    private void processInternally(AminoAcid aminoAcid,
                                   ProteinGraphCalculations fullPlipCalculate,
                                   ProteinGraphCalculations hydrogenPlipCalculate,
                                   ProteinGraphCalculations hydrophobicPlipCalculate,
                                   ProteinGraphCalculations conventionalCalculate) {
        ResidueTopologicProperties fullPlipResidueTopologicProperties = new ResidueTopologicProperties(this,
                fullPlipCalculate.betweenness(aminoAcid),
                fullPlipCalculate.closeness(aminoAcid),
                fullPlipCalculate.clusteringCoefficient(aminoAcid),
                fullPlipCalculate.distinctNeighborhoodCount(aminoAcid));
        ResidueTopologicProperties hydrogenPlipResidueTopologicProperties = new ResidueTopologicProperties(this,
                hydrogenPlipCalculate.betweenness(aminoAcid),
                hydrogenPlipCalculate.closeness(aminoAcid),
                hydrogenPlipCalculate.clusteringCoefficient(aminoAcid),
                hydrogenPlipCalculate.distinctNeighborhoodCount(aminoAcid));
        ResidueTopologicProperties hydrophobicPlipResidueTopologicProperties = new ResidueTopologicProperties(this,
                hydrophobicPlipCalculate.betweenness(aminoAcid),
                hydrophobicPlipCalculate.closeness(aminoAcid),
                hydrophobicPlipCalculate.clusteringCoefficient(aminoAcid),
                hydrophobicPlipCalculate.distinctNeighborhoodCount(aminoAcid));
        ResidueTopologicProperties conventionalResidueTopologicProperties = new ResidueTopologicProperties(this,
                conventionalCalculate.betweenness(aminoAcid),
                conventionalCalculate.closeness(aminoAcid),
                conventionalCalculate.clusteringCoefficient(aminoAcid),
                conventionalCalculate.distinctNeighborhoodCount(aminoAcid));
        aminoAcid.getFeatureContainer().addFeature(new ResidueTopologicPropertiesContainer(this,
                fullPlipResidueTopologicProperties,
                hydrogenPlipResidueTopologicProperties,
                hydrophobicPlipResidueTopologicProperties,
                conventionalResidueTopologicProperties));
    }
}
