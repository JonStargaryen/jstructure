package de.bioforscher.jstructure.si.explorer;

import de.bioforscher.jstructure.si.explorer.model.ExplorerChain;
import de.bioforscher.testutil.TestUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import javax.annotation.PostConstruct;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;

@Service
public class ExplorerService {
    private static final Logger logger = LoggerFactory.getLogger(ExplorerService.class);
    private Map<String, ExplorerChain> chains;
    private List<String> ids;

    @Autowired
    public ExplorerService() {

    }

    @PostConstruct
    public void activate() throws IOException {
        logger.info("starting protein service");

        this.chains = TestUtils.getResourceAsStream("data/pancsa.list")
                .map(line -> line.split(";"))
                .map(ExplorerChain::new)
                .collect(Collectors.toMap(ExplorerChain::getStfId, Function.identity()));

        this.ids = chains.values()
                .stream()
                .map(ExplorerChain::getStfId)
                .sorted()
                .collect(Collectors.toList());
    }

    public ExplorerChain getChain(String stfId) {
        return chains.get(stfId);
    }

    public List<String> getChainIds() {
        return ids;
    }
}