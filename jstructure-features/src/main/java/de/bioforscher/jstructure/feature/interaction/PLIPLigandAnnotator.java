package de.bioforscher.jstructure.feature.interaction;

import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import org.jsoup.nodes.Document;

import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

//TODO infers with other annotator for intra-molecular interactions
public class PLIPLigandAnnotator extends FeatureProvider {
    @Override
    protected void processInternally(Structure protein) {
        protein.chainsWithAminoAcids()
                .parallel()
                .forEach(this::processInternally);

        List<PLIPInteraction> plipInteractions = protein.chainsWithAminoAcids()
                .map(aminoAcid -> aminoAcid.getFeature(PLIPInteractionContainer.class))
                .map(PLIPInteractionContainer::getInteractions)
                .flatMap(Collection::stream)
                .collect(Collectors.toList());

        protein.getFeatureContainer().addFeature(new PLIPInteractionContainer(this, plipInteractions));
    }

    public static Document getDocument(Chain chain) {
        return PLIPRestServiceQuery.calculateLigandDocument(chain);
    }

    private PLIPInteractionContainer processInternally(Chain chain, Document document) {
        return new PLIPInteractionContainer(this, PLIPLigandParser.parse(chain, document));
    }

    private void processInternally(Chain chain) {
        Document document = getDocument(chain);
        PLIPInteractionContainer container = processInternally(chain, document);
        chain.getFeatureContainer().addFeature(container);
        chain.aminoAcids()
                .forEach(aminoAcid -> aminoAcid.getFeatureContainer().addFeature(container.getInteractionsFor(aminoAcid)));
    }
}
