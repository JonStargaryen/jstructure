package de.bioforscher.explorer.membrane.model;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.parser.sifts.SiftsParser;
import de.bioforscher.jstructure.parser.uniprot.UniProtAnnotationContainer;
import de.bioforscher.jstructure.parser.uniprot.UniProtAnnotator;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Represent the multi-sequence alignment of 1 cluster.
 * Created by bittrich on 3/17/17.
 */
public class ExplorerAlignment {
    private static final UniProtAnnotator uniProtAnnotator = new UniProtAnnotator();

    private String representativeChainId;
    private List<String> homologous;
    private List<ExplorerSequencePosition> positions;
    private List<ExplorerSequence> chains;
    private int length;

    public ExplorerAlignment() {
    }

    public ExplorerAlignment(String representativeChainId, List<Chain> chains, Map<String, String> alignedSequences) {
        this.representativeChainId = representativeChainId;
        this.homologous = alignedSequences.keySet().stream().collect(Collectors.toList());
        this.length = alignedSequences.get(representativeChainId).length();

        this.chains = chains.stream()
                .map(ExplorerSequence::new)
                .collect(Collectors.toList());

        // gather uniprot data
        Set<String> uniprotIds = chains.stream()
                .map(chain -> chain.getFeature(String.class, SiftsParser.UNIPROT_ID))
                .collect(Collectors.toSet());

        Set<UniProtAnnotationContainer> uniProtAnnotationContainers = uniprotIds.stream()
                .map(uniProtAnnotator::process)
                .collect(Collectors.toSet());

        this.positions = IntStream.range(0, length)
                .mapToObj(position -> new ExplorerSequencePosition(uniProtAnnotationContainers, alignedSequences, position))
                .collect(Collectors.toList());
    }

    public String getRepresentativeChainId() {
        return representativeChainId;
    }

    public List<ExplorerSequencePosition> getPositions() {
        return positions;
    }

    public List<String> getHomologous() {
        return homologous;
    }

    public int getLength() {
        return length;
    }

    public List<ExplorerSequence> getChains() {
        return chains;
    }
}
