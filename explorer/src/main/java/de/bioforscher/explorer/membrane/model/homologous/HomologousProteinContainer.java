package de.bioforscher.explorer.membrane.model.homologous;

import de.bioforscher.jstructure.parser.ProteinParser;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Information on homologous proteins for a given reference protein.
 * Created by bittrich on 3/15/17.
 */
public class HomologousProteinContainer {
    private static final Logger logger = LoggerFactory.getLogger(HomologousProteinContainer.class);
    private List<String> ids;
    private List<HomologousProtein> proteins;

    public HomologousProteinContainer() {
    }

    public HomologousProteinContainer(List<String> ids) {
        this.ids = ids;
        this.proteins = ids.stream()
                // we parse protein data, so all sequence
                .map(id -> {
                    logger.info("fetching {}", id);
                    return ProteinParser.source(id).parse();
                })
                .map(HomologousProtein::new)
                .collect(Collectors.toList());
    }

    public List<String> getIds() {
        return ids;
    }

    public List<HomologousProtein> getProteins() {
        return proteins;
    }
}
