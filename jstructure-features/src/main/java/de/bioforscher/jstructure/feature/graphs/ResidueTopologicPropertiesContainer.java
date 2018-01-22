package de.bioforscher.jstructure.feature.graphs;

import de.bioforscher.jstructure.model.feature.DefaultFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.feature.FeatureProvider;

@DefaultFeatureProvider(TopologicPropertyCalculator.class)
public class ResidueTopologicPropertiesContainer extends FeatureContainerEntry {
    private final ResidueTopologicProperties fullPlip;
    private final ResidueTopologicProperties hydrogenPlip;
    private final ResidueTopologicProperties hydrophobicPlip;
    private final ResidueTopologicProperties conventional;

    public ResidueTopologicPropertiesContainer(FeatureProvider featureProvider,
                                               ResidueTopologicProperties fullPlip,
                                               ResidueTopologicProperties hydrogenPlip,
                                               ResidueTopologicProperties hydrophobicPlip,
                                               ResidueTopologicProperties conventional) {
        super(featureProvider);
        this.fullPlip = fullPlip;
        this.hydrogenPlip = hydrogenPlip;
        this.hydrophobicPlip = hydrophobicPlip;
        this.conventional = conventional;
    }

    public ResidueTopologicProperties getFullPlip() {
        return fullPlip;
    }

    public ResidueTopologicProperties getHydrogenPlip() {
        return hydrogenPlip;
    }

    public ResidueTopologicProperties getHydrophobicPlip() {
        return hydrophobicPlip;
    }

    public ResidueTopologicProperties getConventional() {
        return conventional;
    }
}
