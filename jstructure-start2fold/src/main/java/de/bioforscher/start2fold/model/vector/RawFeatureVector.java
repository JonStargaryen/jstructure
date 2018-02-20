package de.bioforscher.start2fold.model.vector;

import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceArea;
import de.bioforscher.jstructure.feature.energyprofile.EgorAgreement;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfile;
import de.bioforscher.jstructure.feature.graphs.ResidueTopologicPropertiesContainer;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.util.stream.Collectors;

public class RawFeatureVector extends FeatureVector {
    public RawFeatureVector(double secondaryStructureElementSize, double localHydrogen, double localHydrophobic, double localBackbone, double localInteractions, double nonLocalHydrogen, double nonLocalHydrophobic, double nonLocalBackbone, double nonLocalInteractions, double energy, double egor, double rasa, double betweenness, double closeness, double clusteringCoefficient, double hydrogenBetweenness, double hydrogenCloseness, double hydrogenClusteringCoefficient, double hydrophobicBetweenness, double hydrophobicCloseness, double hydrophobicClusteringCoefficient, double convBetweenness, double convCloseness, double convClusteringCoefficient, double distinctNeighborhoods, double convDistinctNeighborhoods) {
        super(secondaryStructureElementSize, localHydrogen, localHydrophobic, localBackbone, localInteractions, nonLocalHydrogen, nonLocalHydrophobic, nonLocalBackbone, nonLocalInteractions, energy, egor, rasa, betweenness, closeness, clusteringCoefficient, hydrogenBetweenness, hydrogenCloseness, hydrogenClusteringCoefficient, hydrophobicBetweenness, hydrophobicCloseness, hydrophobicClusteringCoefficient, convBetweenness, convCloseness, convClusteringCoefficient, distinctNeighborhoods, convDistinctNeighborhoods);
    }

    public static void assignRawFeatureVector(AminoAcid aminoAcid) {
        GenericSecondaryStructure sse = aminoAcid.getFeature(GenericSecondaryStructure.class);

        PLIPInteractionContainer plipInteractionContainer = aminoAcid.getFeature(PLIPInteractionContainer.class);
        PLIPInteractionContainer nonLocalPlipInteractionContainer = new PLIPInteractionContainer(null,
                plipInteractionContainer
                        .getInteractions()
                        .stream()
                        // interactions have to be non-local
                        .filter(inter -> Math.abs(inter.getPartner1().getResidueIndex() - inter.getPartner2().getResidueIndex()) > 5)
                        .collect(Collectors.toList()));
        PLIPInteractionContainer localPlipInteractionContainer = new PLIPInteractionContainer(null,
                plipInteractionContainer
                        .getInteractions()
                        .stream()
                        // interactions have to be local
                        .filter(inter -> !nonLocalPlipInteractionContainer.getInteractions().contains(inter))
                        .collect(Collectors.toList()));

        ResidueTopologicPropertiesContainer residueTopologicPropertiesContainer =
                aminoAcid.getFeature(ResidueTopologicPropertiesContainer.class);

        // assign features to smooth
        RawFeatureVector featureVector = new RawFeatureVector(sse.getSurroundingSecondaryStructureElement(aminoAcid).getSize(),
                localPlipInteractionContainer.getHydrogenBonds().size(),
                localPlipInteractionContainer.getHydrophobicInteractions().size(),
                localPlipInteractionContainer.getBackboneInteractions().size(),
                localPlipInteractionContainer.getInteractions().size(),

                nonLocalPlipInteractionContainer.getHydrogenBonds().size(),
                nonLocalPlipInteractionContainer.getHydrophobicInteractions().size(),
                nonLocalPlipInteractionContainer.getBackboneInteractions().size(),
                nonLocalPlipInteractionContainer.getInteractions().size(),

                aminoAcid.getFeature(EnergyProfile.class).getSolvationEnergy(),
                aminoAcid.getFeature(EgorAgreement.class).getEgorPrediction(),

                aminoAcid.getFeature(AccessibleSurfaceArea.class).getRelativeAccessibleSurfaceArea(),

                residueTopologicPropertiesContainer.getFullPlip().getBetweenness(),
                residueTopologicPropertiesContainer.getFullPlip().getCloseness(),
                residueTopologicPropertiesContainer.getFullPlip().getClusteringCoefficient(),
                residueTopologicPropertiesContainer.getHydrogenPlip().getBetweenness(),
                residueTopologicPropertiesContainer.getHydrogenPlip().getCloseness(),
                residueTopologicPropertiesContainer.getHydrogenPlip().getClusteringCoefficient(),
                residueTopologicPropertiesContainer.getHydrophobicPlip().getBetweenness(),
                residueTopologicPropertiesContainer.getHydrophobicPlip().getCloseness(),
                residueTopologicPropertiesContainer.getHydrophobicPlip().getClusteringCoefficient(),
                residueTopologicPropertiesContainer.getConventional().getBetweenness(),
                residueTopologicPropertiesContainer.getConventional().getCloseness(),
                residueTopologicPropertiesContainer.getConventional().getClusteringCoefficient(),
                residueTopologicPropertiesContainer.getFullPlip().getDistinctNeighborhoodCount(),
                residueTopologicPropertiesContainer.getConventional().getDistinctNeighborhoodCount());

        aminoAcid.getFeatureContainer().addFeature(featureVector);
    }
}
