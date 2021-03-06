package de.bioforscher.jstructure.feature.interaction;

import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import org.jsoup.nodes.Document;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Annotates residue-residue interactions within a protein.
 * Created by bittrich on 2/9/17.
 */
public class PLIPIntraMolecularAnnotator extends FeatureProvider {
    private static final Logger logger = LoggerFactory.getLogger(PLIPIntraMolecularAnnotator.class);

    @Override
    protected void processInternally(Structure protein) {
        protein.chainsWithAminoAcids()
//                .parallel()
                .forEach(chain -> process(chain, getDocument(chain)));

        List<PLIPInteraction> plipInteractions = protein.chainsWithAminoAcids()
                .map(chain -> chain.getFeature(PLIPInteractionContainer.class))
                .map(PLIPInteractionContainer::getInteractions)
                .flatMap(Collection::stream)
                .collect(Collectors.toList());

        protein.getFeatureContainer().addFeature(new PLIPInteractionContainer(this, plipInteractions));
    }

    public void process(Chain chain) {
        process(chain, getDocument(chain));

        List<PLIPInteraction> plipInteractions = chain.getFeature(PLIPInteractionContainer.class).getInteractions();

        chain.getParentStructure().getFeatureContainer().addFeature(new PLIPInteractionContainer(this, plipInteractions));
    }

    public static Document getDocument(Chain chain) {
        return PLIPRestServiceQuery.getIntraMolecularDocument(chain);
    }

    private PLIPInteractionContainer processInternally(Chain chain, Document document) {
        logger.debug("processing chain '{}_{}' by PLIP",
                chain.getParentStructure().getProteinIdentifier(),
                chain.getChainIdentifier().getChainId());
        return new PLIPInteractionContainer(this, PLIPParser.parse(chain, document));
    }

    public void process(Chain chain, Document document) {
        logger.debug("assigning results from document to {}",
                chain.getChainIdentifier().getFullName());
        PLIPInteractionContainer container = processInternally(chain, document);
        chain.getFeatureContainer().addFeature(container);
        chain.aminoAcids()
                .forEach(aminoAcid -> aminoAcid.getFeatureContainer().addFeature(container.getInteractionsFor(aminoAcid)));
    }
}
