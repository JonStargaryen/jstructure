package de.bioforscher.jstructure.si.explorer;

import de.bioforscher.testutil.TestUtils;

import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class DataSource {
    private static final DataSource INSTANCE = new DataSource();
    private final Map<String, ExplorerChain> chains;

    private DataSource() {
        this.chains = TestUtils.getResourceAsStream("data/ids.list")
                .map(line -> line.split(";"))
//                .filter(split -> split[0].equals("1tup_A"))
                .map(ExplorerChain::new)
                .collect(Collectors.toMap(ExplorerChain::getStfId, Function.identity()));
    }

    public static DataSource getInstance() {
        return INSTANCE;
    }

    public Map<String, ExplorerChain> getChainMap() {
        return chains;
    }

    public Stream<ExplorerChain> chains() {
        return chains.values().stream();
    }

    public Stream<ExplorerChain> start2FoldChains() {
        return chains.values()
                .stream()
                .filter(explorerChain -> explorerChain.getStfId().startsWith("STF"));
    }
}
