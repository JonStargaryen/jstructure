package de.bioforscher.jstructure.feature.interactions;

import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import org.jsoup.nodes.Document;

import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Annotates residue-residue interactions within a protein.
 * Created by bittrich on 2/9/17.
 */
public class PLIPIntraMolecularAnnotator extends FeatureProvider {
    @Override
    protected void processInternally(Structure protein) {
        protein.chainsWithAminoAcids()
                .parallel()
                .forEach(chain -> process(chain, getDocument(chain)));

        List<PLIPInteraction> plipInteractions = protein.chainsWithAminoAcids()
                .map(chain -> chain.getFeature(PLIPInteractionContainer.class))
                .map(PLIPInteractionContainer::getInteractions)
                .flatMap(Collection::stream)
                .collect(Collectors.toList());

        protein.getFeatureContainer().addFeature(new PLIPInteractionContainer(this, plipInteractions));
    }

    public static Document getDocument(Chain chain) {
        return PLIPRestServiceQuery.getIntraMolecularDocument(chain);
    }

    private PLIPInteractionContainer processInternally(Chain chain, Document document) {
        return new PLIPInteractionContainer(this, PLIPParser.parse(chain, document));
    }

    public void process(Chain chain, Document document) {
        PLIPInteractionContainer container = processInternally(chain, document);
        chain.getFeatureContainer().addFeature(container);
        chain.aminoAcids()
                .forEach(aminoAcid -> aminoAcid.getFeatureContainer().addFeature(container.getInteractionsFor(aminoAcid)));
    }
}
