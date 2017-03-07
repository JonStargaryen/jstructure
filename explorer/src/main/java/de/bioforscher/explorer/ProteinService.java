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
import de.bioforscher.jstructure.parser.kinkfinder.KinkFinderParser;
import de.bioforscher.jstructure.parser.plip.PLIPAnnotator;
import de.bioforscher.jstructure.parser.uniprot.UniProtAnnotator;
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
import java.util.concurrent.Callable;
import java.util.concurrent.Executors;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Access to the database service.
 * Created by bittrich on 2/20/17.
 */
@Service
public class ProteinService {
    private static final Logger logger = LoggerFactory.getLogger(ProteinService.class);
    private ProteinRepository repository;
    private List<String> allProteins;
    private List<String> nonRedundantAlphaHelicalProteinIds;
    private List<String> features = Stream.of(UniProtAnnotator.UNIPROT_ANNOTATION,
            AccessibleSurfaceAreaCalculator.RELATIVE_ACCESSIBLE_SURFACE_AREA,
            SequenceMotifAnnotator.SEQUENCE_MOTIF,
            SecondaryStructureAnnotator.SECONDARY_STRUCTURE_STATES,
            ANVIL.MEMBRANE,
            PLIPAnnotator.PLIP_INTERACTIONS).collect(Collectors.toList());

    @Autowired
    public ProteinService(ProteinRepository repository) {
        this.repository = repository;
    }

    @PostConstruct
    public void activate() throws IOException {
        // clear database
        repository.deleteAll();

        // load protein lists
        List<String> allProteinIds = Files.lines(getResource("pdbtm_all.list")) //http://pdbtm.enzim.hu/data/pdbtm_all.list
                .distinct()
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
//        List<String> proteinsToProcess = this.nonRedundantAlphaHelicalProteinIds.stream()
//                .map(line -> line.substring(0, 4))
//                .distinct()
//                .limit(0)
//                .filter(id -> !knownProteins.contains(id))
//                .collect(Collectors.toList());

        List<String> proteinsToProcess = Stream.of("5A63", "1BRR", "4NEF").collect(Collectors.toList());
        this.allProteins = proteinsToProcess;

        proteinsToProcess.stream()
                .map(pdbid -> (Callable) () -> {
                    try{
                        logger.info("fetching information for {}", pdbid);
                        Protein protein = ProteinParser.parseProteinById(pdbid);
                        featureProviders.forEach(abstractFeatureProvider -> abstractFeatureProvider.process(protein));
                        KinkFinderParser.parseKinkFinderFile(protein, Paths.get("/home/bittrich/git/phd_sb_repo/data/kink_finder/" + protein.getName().toLowerCase() + ".kinks"));
                        repository.save(new ExplorerProtein(protein));
                    } catch (Exception e) {
                        logger.error("gathering information for {} failed: {}", pdbid, e.getLocalizedMessage());
                    }
                    return null;
                })
                .forEach(task -> Executors.newWorkStealingPool().submit(task));

//        this.allProteins = repository.findAll().stream()
//                .map(ExplorerProtein::getName)
//                .collect(Collectors.toList());
    }

    public List<String> getAllProteins() {
        return allProteins;
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
