package de.bioforscher.explorer.model;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.parser.sifts.ChainSiftsMapping;
import de.bioforscher.jstructure.parser.sifts.SiftsMappingProvider;
import de.bioforscher.jstructure.parser.uniprot.UniProtAnnotationContainer;
import de.bioforscher.jstructure.parser.uniprot.UniProtAnnotator;

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
                .map(chain -> chain.getParentProtein().getName().toLowerCase() + "_" + chain.getChainId())
                .collect(Collectors.toList());

        this.chains = originalChains.stream()
                .map(ExplorerSequence::new)
                .collect(Collectors.toList());

        Set<String> uniProtIds = originalChains.stream()
                .map(chain -> chain.getFeature(ChainSiftsMapping.class, SiftsMappingProvider.SIFTS_MAPPING))
                .map(ChainSiftsMapping::getUniProtId)
                .collect(Collectors.toSet());

        Set<UniProtAnnotationContainer> containers = uniProtIds.stream()
                .map(uniProtAnnotator::process)
                .collect(Collectors.toSet());

        this.length = containers.stream()
                .map(UniProtAnnotationContainer::getSequence)
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
