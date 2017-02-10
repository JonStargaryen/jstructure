package de.bioforscher.jstructure.parser.plip;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.UncheckedIOException;
import java.util.ArrayList;
import java.util.List;

import static de.bioforscher.jstructure.parser.plip.PLIPAnnotator.PLIP_INTERACTIONS;

/**
 * Annotates residue-residue interactions within a protein.
 * Created by bittrich on 2/9/17.
 */
@FeatureProvider(provides = PLIP_INTERACTIONS)
public class PLIPAnnotator extends AbstractFeatureProvider {
    private static final Logger logger = LoggerFactory.getLogger(PLIPAnnotator.class);
    public static final String PLIP_INTERACTIONS = "PLIP_INTERACTIONS";

    @Override
    protected void processInternally(Protein protein) {
        List<PLIPInteraction> globalPlipInteractions = new ArrayList<>();

        protein.chains()
                // for safety: ignore non-amino-acid chains
                .filter(chain -> chain.aminoAcids().count() > 0)
                .forEach((Chain chain) -> {
                        try {
                            String plipXmlContent = PLIPRestServiceQuery.getPlipResults(protein.getName(), chain.getChainId());
                            PLIPParser.parse(chain, plipXmlContent).forEach(plipInteraction -> {
                                Group group = plipInteraction.getPartner1();
                                assignInteraction(group, plipInteraction);
                                checkedAddInteraction(globalPlipInteractions, plipInteraction);
                            });
                        } catch (UncheckedIOException e) {
                            logger.warn("failed to fetch PLIP results");
                        }
        });

        protein.setFeature(PLIP_INTERACTIONS, globalPlipInteractions);
    }

    private void checkedAddInteraction(List<PLIPInteraction> globalPlipInteractions, PLIPInteraction plipInteractionToAdd) {
        // do not store symmetric interactions, but only one direction for the global interaction annotation
        PLIPInteractionType plipInteractionType = plipInteractionToAdd.getPlipInteractionType();
        Group partner1 = plipInteractionToAdd.getPartner1();
        Group partner2 = plipInteractionToAdd.getPartner2();
        if(globalPlipInteractions.stream()
                .filter(plipInteraction -> plipInteraction.getPlipInteractionType().equals(plipInteractionType))
                .anyMatch(plipInteraction -> plipInteraction.getPartner1().equals(partner2) && plipInteraction.getPartner2().equals(partner1))) {
            return;
        }

        globalPlipInteractions.add(plipInteractionToAdd);
    }

    @SuppressWarnings("unchecked")
    private void assignInteraction(Group group, PLIPInteraction plipInteraction) {
        List<PLIPInteraction> value = group.getFeature(List.class, PLIP_INTERACTIONS);
        // entry will be null at first - create list and assign reference
        if(value == null) {
            value = new ArrayList<>();
            group.setFeature(PLIP_INTERACTIONS, value);
        }
        value.add(plipInteraction);
    }
}
