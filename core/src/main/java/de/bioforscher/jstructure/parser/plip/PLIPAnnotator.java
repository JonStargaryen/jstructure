package de.bioforscher.jstructure.parser.plip;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.plip.interaction.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.UncheckedIOException;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import static de.bioforscher.jstructure.parser.plip.PLIPAnnotator.PLIP_INTERACTIONS;

/**
 * Annotates residue-residue interactions within a protein.
 * TODO this implementation is really ugly :x
 * Created by bittrich on 2/9/17.
 */
@FeatureProvider(provides = { PLIP_INTERACTIONS })
public class PLIPAnnotator extends AbstractFeatureProvider {
    private static final Logger logger = LoggerFactory.getLogger(PLIPAnnotator.class);
    public static final String PLIP_INTERACTIONS = "PLIP_INTERACTIONS";

    @Override
    protected void processInternally(Protein protein) {
        List<PLIPInteraction> globalPlipInteractions = new ArrayList<>();

        protein.chains()
                // for safety: ignore non-amino-acid chains
                .filter(chain -> chain.aminoAcids().count() > 0)
                .forEach(chain -> {
                    try {
                        String plipXmlContent = PLIPRestServiceQuery.getPlipResults(protein.getName(), chain.getChainId());
                        List<PLIPInteraction> plipInteractions = PLIPParser.parse(chain, plipXmlContent);
                        globalPlipInteractions.addAll(plipInteractions);
                    } catch (UncheckedIOException e) {
                        logger.warn("failed to fetch PLIP results");
                    }
        });

        List<HalogenBond> halogenBonds = filterInteractions(globalPlipInteractions, HalogenBond.class);
        List<HydrogenBond> hydrogenBonds = filterInteractions(globalPlipInteractions, HydrogenBond.class);
        List<HydrophobicInteraction> hydrophobicInteractions = filterInteractions(globalPlipInteractions, HydrophobicInteraction.class);
        List<MetalComplex> metalComplexes = filterInteractions(globalPlipInteractions, MetalComplex.class);
        List<PiCationInteraction> piCationInteractions = filterInteractions(globalPlipInteractions, PiCationInteraction.class);
        List<PiStacking> piStackings = filterInteractions(globalPlipInteractions, PiStacking.class);
        List<SaltBridge> saltBridges = filterInteractions(globalPlipInteractions, SaltBridge.class);
        List<WaterBridge> waterBridges = filterInteractions(globalPlipInteractions, WaterBridge.class);

        protein.setFeature(PLIP_INTERACTIONS, new PLIPInteractionContainer(halogenBonds, hydrogenBonds,
                hydrophobicInteractions, metalComplexes, piCationInteractions, piStackings, saltBridges, waterBridges));
    }

    private <I> List<I> filterInteractions(List<PLIPInteraction> globalPlipInteractions, Class<I> content) {
        return globalPlipInteractions.stream()
                .filter(content::isInstance)
                .map(content::cast)
                .collect(Collectors.toList());
    }
}
