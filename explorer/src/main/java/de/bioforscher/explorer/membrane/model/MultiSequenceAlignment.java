package de.bioforscher.explorer.membrane.model;

import de.bioforscher.jstructure.model.Pair;

import java.util.List;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * Represent the multi-sequence alignemnt of 1 cluster.
 * Created by bittrich on 3/17/17.
 */
public class MultiSequenceAlignment {
    private String representativeChainId;
    private List<ExplorerSequence> chains;

    public MultiSequenceAlignment() {
    }

    public MultiSequenceAlignment(String representativeChainId, List<ExplorerChain> explorerChains, String alignment) {
        this.representativeChainId = representativeChainId;

        List<String> split = Pattern.compile(">")
                .splitAsStream(alignment)
                .skip(1)
                .collect(Collectors.toList());

        this.chains = split.stream()
                .map(s -> s.split("SEQUENCE"))
                // substring for pdbId:chainId
//                .map(s -> new Pair<>(s[0].substring(0, 6), Pattern.compile("[\\r\\n]+")
//                        .splitAsStream(s[1])
//                        .collect(Collectors.joining())))
                .map(s ->  new Pair<>(s[0].substring(0, 6), s[1].replaceAll("\\s", "")))
                .map(pair -> new ExplorerSequence(resolveChain(explorerChains, pair.getLeft()), pair.getLeft(), pair.getRight()))
                .collect(Collectors.toList());
    }

    private ExplorerChain resolveChain(List<ExplorerChain> explorerChains, String id) {
        return explorerChains.stream()
                .filter(chain -> chain.getId().equals(id))
                .findFirst()
                .get();
    }

    public String getRepresentativeChainId() {
        return representativeChainId;
    }

    public List<ExplorerSequence> getChains() {
        return chains;
    }
}
