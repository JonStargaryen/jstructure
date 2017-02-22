package de.bioforscher.explorer;

import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import javax.annotation.PostConstruct;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;

/**
 * Access to the database service.
 * Created by bittrich on 2/20/17.
 */
@Service
public class ProteinService {
    private final ProteinRepository repository;
    private List<String> allProteinIds;
    private List<String> nonRedundantAlphaHelicalProteinIds;

    @Autowired
    public ProteinService(ProteinRepository repository) throws IOException {
        this.repository = repository;
    }

    @PostConstruct
    public void activate() throws IOException {
        repository.deleteAll();

        // load protein lists
        this.allProteinIds = Files.lines(getResource("pdbtm_all.list")) //http://pdbtm.enzim.hu/data/pdbtm_all.list
                .limit(1)
                .map(line -> line.substring(0, 4))
                .map(String::toUpperCase)
                .collect(Collectors.toList());
        this.nonRedundantAlphaHelicalProteinIds = Files.lines(getResource("pdbtm_alpha_nr.list")) //http://pdbtm.enzim.hu/data/pdbtm_alpha_nr.list
                .map(String::toUpperCase)
                .collect(Collectors.toList());

        // parse protein data
        allProteinIds.forEach(pdbid -> {
                    Protein protein = ProteinParser.parseProteinById(pdbid);
                    repository.save(protein);
                });
    }

    public List<String> getAllProteins() {
        return allProteinIds;
    }

    public List<String> getNonRedundantAlphaHelicalProteins() {
        return nonRedundantAlphaHelicalProteinIds;
    }

    public Protein getProtein(String pdbid) {
        return repository.findByTheProteinsName(pdbid.toUpperCase()).get(0);
    }

    private Path getResource(String filename) {
        // some a bit hacky way to ensure correct paths on windows (as some / will be added as prefix)
        return Paths.get(Objects.requireNonNull(Thread.currentThread().getContextClassLoader().getResource(filename)).getPath().replaceFirst("^/(.:/)", "$1"));
    }
}
