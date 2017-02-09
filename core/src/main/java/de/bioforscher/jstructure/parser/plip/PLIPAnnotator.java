package de.bioforscher.jstructure.parser.plip;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;

import java.util.ArrayList;
import java.util.List;

import static de.bioforscher.jstructure.parser.plip.PLIPAnnotator.PLIP_INTERACTIONS;

/**
 * Annotates residue-residue interactions within a protein.
 * Created by bittrich on 2/9/17.
 */
@FeatureProvider(provides = PLIP_INTERACTIONS)
public class PLIPAnnotator extends AbstractFeatureProvider {
    public static final String PLIP_INTERACTIONS = "PLIP_INTERACTIONS";

    @Override
    protected void processInternally(Protein protein) {
        protein.chains()
                // for safety: ignore non-amino-acid chains
                .filter(chain -> chain.aminoAcids().count() > 0)
                .forEach((Chain chain) -> {
            String plipXmlContent = PLIPRestServiceQuery.getPlipResults(protein.getName(), chain.getChainId());
            PLIPParser.parse(chain, plipXmlContent).forEach(plipInteraction -> {
                Group group = plipInteraction.getPartner1();
                assignInteraction(group, plipInteraction);
            });
        });
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
