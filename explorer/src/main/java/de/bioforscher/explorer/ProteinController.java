package de.bioforscher.explorer;

import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.web.bind.annotation.PathVariable;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RequestMethod;
import org.springframework.web.bind.annotation.RestController;

import java.io.IOException;
import java.util.List;

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
    public List<String> getAllProteins() throws IOException {
        return proteinService.getAllProteins();
    }

    @RequestMapping(value = "/alpha_nr", method = RequestMethod.GET)
    public List<String> getNonRedundantAlphaHelicalProteins() throws IOException {
        return proteinService.getNonRedundantAlphaHelicalProteins();
    }

    @RequestMapping(value = "/pdbid/{pdbid}", method = RequestMethod.GET)
    public Protein getProteinByID(@PathVariable("pdbid") String pdbid) {
        return ProteinParser.parseProteinById(pdbid);
    }
}
