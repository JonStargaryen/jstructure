package de.bioforscher.jstructure.si.explorer;

import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.web.bind.annotation.PathVariable;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RequestMethod;
import org.springframework.web.bind.annotation.RestController;

import java.io.IOException;
import java.util.List;
import java.util.NoSuchElementException;

@RestController
@RequestMapping(value = "/api/", method = RequestMethod.GET)
public class ExplorerController {
    private final ExplorerService explorerService;

    @Autowired
    public ExplorerController(ExplorerService explorerService) throws IOException {
        this.explorerService = explorerService;
    }

    @RequestMapping(value = "/ids", method = RequestMethod.GET)
    public List<String> getChainIds() {
        return explorerService.getChainIds();
    }

    @RequestMapping(value = "/json/{id}", method = RequestMethod.GET)
    public ExplorerChain getChainById(@PathVariable("id") String id) throws NoSuchElementException {
        return explorerService.getChain(id);
    }

    @RequestMapping(value = "/pdb/{id}", method = RequestMethod.GET, produces = "text/plain")
    public String getStructureById(@PathVariable("id") String id) throws NoSuchElementException {
        return explorerService.getChain(id).getPdbRepresentation();
    }
}
