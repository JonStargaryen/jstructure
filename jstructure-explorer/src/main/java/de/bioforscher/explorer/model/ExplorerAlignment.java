package de.bioforscher.explorer.model;

import de.bioforscher.jstructure.feature.mapping.ChainMapping;
import de.bioforscher.jstructure.feature.uniprot.UniProtAnnotationContainer;
import de.bioforscher.jstructure.feature.uniprot.UniProtAnnotator;
import de.bioforscher.jstructure.model.structure.Chain;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * The sequence alignment of one cluster.
 * Created by bittrich on 4/6/17.
 */
@Deprecated
public class ExplorerAlignment {
    private static final UniProtAnnotator uniProtAnnotator = new UniProtAnnotator();
    private String representativeId;
    private List<String> homologous;
    private List<ExplorerSequencePosition> positions;
    private List<ExplorerSequence> chains;
    private int length;

    public ExplorerAlignment() {

    }

    public ExplorerAlignment(List<Chain> originalChains, String representativeId) {
        this.representativeId = representativeId;
        this.homologous = originalChains.stream()
                .map(chain -> chain.getChainId().getFullName())
                .collect(Collectors.toList());

        this.chains = originalChains.stream()
                .map(ExplorerSequence::new)
                .collect(Collectors.toList());

        Set<String> uniProtIds = originalChains.stream()
                .map(chain -> chain.getFeatureContainer().getFeature(ChainMapping.class).getUniProtId())
                .collect(Collectors.toSet());

        Set<UniProtAnnotationContainer> containers = uniProtIds.stream()
                .map(uniProtAnnotator::processUniProtId)
                .collect(Collectors.toSet());

        this.length = containers.stream()
                .map(UniProtAnnotationContainer::getUniProtSequence)
                .mapToInt(String::length)
                .max()
                .orElse(0);

        this.positions = IntStream.range(0, length)
                .mapToObj(position -> new ExplorerSequencePosition(containers, originalChains, position))
                .collect(Collectors.toList());
    }

    public String getRepresentativeId() {
        return representativeId;
    }

    public List<String> getHomologous() {
        return homologous;
    }

    public List<ExplorerSequencePosition> getPositions() {
        return positions;
    }

    public List<ExplorerSequence> getChains() {
        return chains;
    }

    public int getLength() {
        return length;
    }
}
