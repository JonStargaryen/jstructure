package de.bioforscher.explorer;

import de.bioforscher.explorer.model.ExplorerAlignment;
import de.bioforscher.explorer.model.ExplorerChain;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.web.bind.annotation.PathVariable;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RequestMethod;
import org.springframework.web.bind.annotation.RestController;

import java.io.IOException;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * Provides REST access to the stored protein data.
 * Created by bittrich on 2/20/17.
 */
@RestController
@RequestMapping(value = "/api/chains/", method = RequestMethod.GET)
public class ChainController {
    private static final Logger logger = LoggerFactory.getLogger(ChainController.class);
    private final ChainService chainService;

    @Autowired
    public ChainController(ChainService chainService) throws IOException {
        this.chainService = chainService;
    }

    @RequestMapping(value = "/reps", method = RequestMethod.GET)
    public List<String> getAllRepresentativeChainIds() {
        return chainService.getAllRepresentativeChainIds();
    }

    @RequestMapping(value = "/all", method = RequestMethod.GET)
    public List<String> getAllChainIds() {
        return chainService.getAllChainIds();
    }

    @RequestMapping(value = "/json/{id}", method = RequestMethod.GET)
    public ExplorerChain getChainById(@PathVariable("id") String id) throws NoSuchElementException {
        return chainService.getChain(id);
    }

    @RequestMapping(value = "/alignment/{id}", method = RequestMethod.GET)
    public ExplorerAlignment getAlignmentById(@PathVariable("id") String id) throws NoSuchElementException {
        return chainService.getAlignment(id);
    }

    @RequestMapping(value = "/pdb/{id}", method = RequestMethod.GET, produces = "text/plain")
    public String getStructureById(@PathVariable("id") String id) throws NoSuchElementException {
        return chainService.getChain(id).getPdb();
    }
}
