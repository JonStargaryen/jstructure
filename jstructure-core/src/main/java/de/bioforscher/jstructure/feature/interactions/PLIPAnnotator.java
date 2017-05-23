package de.bioforscher.jstructure.feature.interactions;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Protein;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Annotates residue-residue interactions within a protein.
 * Created by bittrich on 2/9/17.
 */
@FeatureProvider(provides = PLIPInteractionContainer.class)
public class PLIPAnnotator extends AbstractFeatureProvider {
    private static final Logger logger = LoggerFactory.getLogger(PLIPAnnotator.class);

    @Override
    protected void processInternally(Protein protein) {
        List<PLIPInteraction> plipInteractions = protein.chainsWithAminoAcids()
                .parallel()
                .flatMap(chain -> PLIPParser.parse(chain, PLIPRestServiceQuery.getPlipResults(chain.getChainId())).stream())
                .collect(Collectors.toList());

        protein.getFeatureContainer().addFeature(new PLIPInteractionContainer(this, plipInteractions));
    }
}
