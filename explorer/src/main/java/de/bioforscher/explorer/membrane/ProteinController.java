package de.bioforscher.explorer.membrane;

import de.bioforscher.explorer.membrane.model.ExplorerProtein;
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
@RequestMapping(value = "/api/proteins/", method = RequestMethod.GET)
public class ProteinController {
    private final ProteinService proteinService;

    @Autowired
    public ProteinController(ProteinService proteinService) throws IOException {
        this.proteinService = proteinService;
    }

    @RequestMapping(value = "/all", method = RequestMethod.GET)
    public List<String> getAllProteins() {
        return proteinService.getAllProteins();
    }

    @RequestMapping(value = "/pdbid/{pdbid}", method = RequestMethod.GET)
    public ExplorerProtein getProteinByID(@PathVariable("pdbid") String pdbId) throws NoSuchElementException {
        return proteinService.getProtein(pdbId);
    }
}
