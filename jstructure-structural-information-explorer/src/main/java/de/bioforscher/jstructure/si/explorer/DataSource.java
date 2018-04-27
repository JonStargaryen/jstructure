package de.bioforscher.jstructure.si.explorer;

import de.bioforscher.testutil.TestUtils;

import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;

public class DataSource {
    private static final DataSource INSTANCE = new DataSource();
    private final Map<String, ExplorerChain> chains;

    private DataSource() {
        this.chains = TestUtils.getResourceAsStream("data/ids.list")
                .map(line -> line.split(";"))
                .map(ExplorerChain::new)
                .collect(Collectors.toMap(ExplorerChain::getStfId, Function.identity()));
    }

    public static DataSource getInstance() {
        return INSTANCE;
    }

    public Map<String, ExplorerChain> getChains() {
        return chains;
    }
}
