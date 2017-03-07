package de.bioforscher.helices;

import de.bioforscher.explorer.model.ExplorerProtein;
import de.bioforscher.helices.model.AnnotatedHelix;
import de.bioforscher.helices.model.HelixClassification;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import javax.annotation.PostConstruct;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Objects;
import java.util.concurrent.Executors;
import java.util.stream.Collectors;

/**
 * Access to the database service.
 * Created by bittrich on 3/7/17.
 */
@Service
public class HelixService {
    private static final Logger logger = LoggerFactory.getLogger(HelixService.class);
    private HelixRepository repository;
    private List<String> helixLines;
    private List<ExplorerProtein> proteins;

    @Autowired
    public HelixService(HelixRepository repository) {
        this.repository = repository;
    }

    @PostConstruct
    public void activate() throws IOException {
        // clear database
        repository.deleteAll();

        // load protein lists
        helixLines = Files.lines(getResource("ahah_gold_standard.csv")) //http://opig.stats.ox.ac.uk/webapps/ahah/php/experiment_results.php
                .filter(line -> !line.startsWith(">") && !line.startsWith("helix_id"))
                .collect(Collectors.toList());

        proteins = helixLines.stream()
                .map(line -> line.substring(0, 4))
                .distinct()
                .map(ExplorerProtein::createForHelixExplorer)
                .collect(Collectors.toList());

        helixLines.forEach(helixLine -> Executors.newWorkStealingPool().submit(() -> {
            String[] helixLineSplit = helixLine.split(",");
            String helixId = helixLineSplit[0];
            String[] idSplit = helixId.split("_");

            logger.info("composing information for {}", helixId);
            ExplorerProtein protein = proteins.stream()
                    .filter(candidate -> candidate.getName().equalsIgnoreCase(helixLine.split("_")[0]))
                    .findFirst()
                    .get();

            AnnotatedHelix annotatedHelix = new AnnotatedHelix(helixId,
                    idSplit[0],
                    idSplit[1],
                    Integer.valueOf(idSplit[2]),
                    Integer.valueOf(idSplit[3]),
                    Double.valueOf(helixLineSplit[2]),
                    Double.valueOf(helixLineSplit[3]),
                    Double.valueOf(helixLineSplit[4]),
                    Integer.valueOf(helixLineSplit[5]),
                    mapToClassification(helixLineSplit[1]),
                    protein);

            repository.save(annotatedHelix);
        }));
    }

    private HelixClassification mapToClassification(String string) {
        switch(string) {
            case "KINK":
                return HelixClassification.KINK;
            case "CURVED":
                return HelixClassification.CURVED;
            case "STRAIGHT":
                return HelixClassification.STRAIGHT;
            case "UNCLASSIFIED":
                return HelixClassification.UNCLASSIFIED;
            default:
                throw new UnsupportedOperationException("do not know how to handle: " + string);
        }
    }

    public List<String> getAllHelices() {
        return helixLines.stream()
                .map(line -> line.split(",")[0])
                .collect(Collectors.toList());
    }

    public AnnotatedHelix getHelix(String helixid) {
        List<AnnotatedHelix> result = repository.findByTheHelixId(helixid.toUpperCase());
        if(result.size() > 0) {
            return result.get(0);
        }
        throw new NoSuchElementException("could not retrieve entry for " + helixid);
    }

    private Path getResource(String filename) {
        // some a bit hacky way to ensure correct paths on windows (as some / will be added as prefix)
        return Paths.get(Objects.requireNonNull(Thread.currentThread().getContextClassLoader().getResource(filename)).getPath().replaceFirst("^/(.:/)", "$1"));
    }
}
