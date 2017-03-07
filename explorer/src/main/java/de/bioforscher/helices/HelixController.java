package de.bioforscher.helices;

import de.bioforscher.helices.model.AnnotatedHelix;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.web.bind.annotation.PathVariable;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RequestMethod;
import org.springframework.web.bind.annotation.RestController;

import java.io.IOException;
import java.util.List;

/**
 * Provides REST access to the stored AHAH data.
 * Created by bittrich on 3/7/17.
 */
@RestController
@RequestMapping(value = "/api/proteins/", method = RequestMethod.GET)
public class HelixController {
    private final HelixService helixService;

    @Autowired
    public HelixController(HelixService helixService) throws IOException {
        this.helixService = helixService;
    }

    @RequestMapping(value = "/all", method = RequestMethod.GET)
    public List<String> getAllHelices() throws IOException {
        return helixService.getAllHelices();
    }

    @RequestMapping(value = "/{helixid}", method = RequestMethod.GET)
    public AnnotatedHelix getHelixByID(@PathVariable("helixid") String helixid) {
        return helixService.getHelix(helixid);
    }
}