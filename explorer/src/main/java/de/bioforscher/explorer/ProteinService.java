package de.bioforscher.explorer;

import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

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
    private final List<String> allProteinIds;
    private final List<String> nonRedundantAlphaHelicalProteinIds;

    @Autowired
    public ProteinService(ProteinRepository repository) throws IOException {
        this.repository = repository;

        // load protein lists
        this.allProteinIds = Files.lines(getResource("pdbtm_all.list")) //http://pdbtm.enzim.hu/data/pdbtm_all.list
                .limit(5)
                .collect(Collectors.toList());
        this.nonRedundantAlphaHelicalProteinIds = Files.lines(getResource("pdbtm_alpha_nr.list")) //http://pdbtm.enzim.hu/data/pdbtm_alpha_nr.list
                .limit(5)
                .collect(Collectors.toList());
    }

    public List<String> getAllProteins() {
        return allProteinIds;
    }

    public List<String> getNonRedundantAlphaHelicalProteins() {
        return nonRedundantAlphaHelicalProteinIds;
    }

    private static Path getResource(String filename) {
        // some a bit hacky way to ensure correct paths on windows (as some / will be added as prefix)
        return Paths.get(Objects.requireNonNull(Thread.currentThread().getContextClassLoader().getResource(filename)).getPath().replaceFirst("^/(.:/)", "$1"));
    }
}
