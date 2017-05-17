package de.bioforscher.jstructure.feature.interactions;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Protein;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.UncheckedIOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Annotates residue-residue interactions within a protein.
 * Created by bittrich on 2/9/17.
 */
@FeatureProvider(provides = PLIPInteractionContainer.class)
public class PLIPAnnotator extends AbstractFeatureProvider {
    private static final Logger logger = LoggerFactory.getLogger(PLIPAnnotator.class);

    @Override
    protected void processInternally(Protein protein) {
        List<PLIPInteraction> globalPlipInteractions = new ArrayList<>();

        protein.aminoAcidChains()
                .forEach(chain -> {
                    try {
                        String plipXmlContent = PLIPRestServiceQuery.getPlipResults(protein.getName(), chain.getChainId());
                        globalPlipInteractions.addAll(PLIPParser.parse(chain, plipXmlContent));
                    } catch (UncheckedIOException e) {
                        logger.warn("failed to fetch PLIP results");
                    }
        });

        protein.getFeatureContainer().addFeature(new PLIPInteractionContainer(this, globalPlipInteractions));
    }
}
