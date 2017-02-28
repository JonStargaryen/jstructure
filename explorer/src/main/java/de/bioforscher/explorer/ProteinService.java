package de.bioforscher.explorer;

import de.bioforscher.explorer.model.ExplorerProtein;
import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceAreaCalculator;
import de.bioforscher.jstructure.feature.motif.SequenceMotifAnnotator;
import de.bioforscher.jstructure.feature.sse.SecondaryStructureAnnotator;
import de.bioforscher.jstructure.feature.topology.ANVIL;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProviderRegistry;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import de.bioforscher.jstructure.parser.plip.PLIPAnnotator;
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
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Access to the database service.
 * Created by bittrich on 2/20/17.
 */
@Service
public class ProteinService {
    private static final Logger logger = LoggerFactory.getLogger(ProteinService.class);
    private final ProteinRepository repository;
    private List<String> nonRedundantAlphaHelicalProteinIds;
    private List<String> features = Stream.of(AccessibleSurfaceAreaCalculator.RELATIVE_ACCESSIBLE_SURFACE_AREA,
            SequenceMotifAnnotator.SEQUENCE_MOTIF,
            ANVIL.MEMBRANE,
            SecondaryStructureAnnotator.SECONDARY_STRUCTURE_STATES,
            PLIPAnnotator.PLIP_INTERACTIONS).collect(Collectors.toList());

    @Autowired
    public ProteinService(ProteinRepository repository) throws IOException {
        this.repository = repository;
    }

    @PostConstruct
    public void activate() throws IOException {
        // clear database
//        repository.deleteAll();

        // load protein lists
        List<String> allProteinIds = Files.lines(getResource("pdbtm_all.list")) //http://pdbtm.enzim.hu/data/pdbtm_all.list
                .distinct()
                .limit(1)
                .map(line -> line.substring(0, 4))
                .map(String::toUpperCase)
                .collect(Collectors.toList());
        this.nonRedundantAlphaHelicalProteinIds = Files.lines(getResource("pdbtm_alpha_nr.list")) //http://pdbtm.enzim.hu/data/pdbtm_alpha_nr.list
                .map(String::toUpperCase)
                .collect(Collectors.toList());

        // resolve all feature providers
        List<AbstractFeatureProvider> featureProviders = features.stream()
                .map(FeatureProviderRegistry::resolve)
                .collect(Collectors.toList());

        List<String> knownProteins = repository.findAll().stream()
                .map(ExplorerProtein::getName)
                .collect(Collectors.toList());

        // parse protein data
        List<String> proteinsToProcess = allProteinIds.stream()
                .filter(id -> !knownProteins.contains(id))
                .collect(Collectors.toList());

//        proteinsToProcess.stream()
//                .map(pdbid -> (Callable) () -> {
//                    Protein protein = ProteinParser.parseProteinById(pdbid);
//                    featureProviders.forEach(abstractFeatureProvider ->abstractFeatureProvider.process(protein));
//                    repository.save(new ExplorerProtein(protein));
//                    logger.info("gathered all information for entry {}",pdbid);
//                    return null;
//                })
//                .forEach(task -> Executors.newFixedThreadPool(1).submit(task));

//        List<String> proteinsToProcess = Stream.of("1a0s", "1brr").collect(Collectors.toList());

        proteinsToProcess.forEach(pdbid -> {
            Protein protein = ProteinParser.parseProteinById(pdbid);
            featureProviders.forEach(abstractFeatureProvider ->abstractFeatureProvider.process(protein));
            repository.save(new ExplorerProtein(protein));
        });
    }

    public List<String> getAllProteins() {
        return repository.findAll().stream()
                .map(ExplorerProtein::getName)
                .collect(Collectors.toList());
    }

    public List<String> getNonRedundantAlphaHelicalProteins() {
        return nonRedundantAlphaHelicalProteinIds;
    }

    public ExplorerProtein getProtein(String pdbid) {
        List<ExplorerProtein> result = repository.findByTheProteinsName(pdbid.toUpperCase());
        if(result.size() > 0) {
            return result.get(0);
        }

        throw new NoSuchElementException("could not retrieve entry for " + pdbid);
    }

    private Path getResource(String filename) {
        // some a bit hacky way to ensure correct paths on windows (as some / will be added as prefix)
        return Paths.get(Objects.requireNonNull(Thread.currentThread().getContextClassLoader().getResource(filename)).getPath().replaceFirst("^/(.:/)", "$1"));
    }
}
