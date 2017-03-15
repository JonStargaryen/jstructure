package de.bioforscher.explorer.membrane;

import de.bioforscher.explorer.membrane.model.representative.ExplorerProtein;
import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceAreaCalculator;
import de.bioforscher.jstructure.feature.cerosene.SequenceCerosene;
import de.bioforscher.jstructure.feature.motif.SequenceMotifAnnotator;
import de.bioforscher.jstructure.feature.sse.SecondaryStructureAnnotator;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProviderRegistry;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import de.bioforscher.jstructure.parser.opm.OPMDatabaseQuery;
import de.bioforscher.jstructure.parser.sifts.SiftsParser;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import javax.annotation.PostConstruct;
import java.io.IOException;
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
    private ProteinRepository repository;
    private List<String> allProteins;
    private List<AbstractFeatureProvider> featureProviders = Stream.of(AccessibleSurfaceAreaCalculator.RELATIVE_ACCESSIBLE_SURFACE_AREA,
            SecondaryStructureAnnotator.SECONDARY_STRUCTURE_STATES,
            SequenceCerosene.SEQUENCE_CEROSENE_REPRESENTATION,
            SequenceMotifAnnotator.SEQUENCE_MOTIF,
            SiftsParser.UNIPROT_ID)
            .map(FeatureProviderRegistry::resolve)
            .collect(Collectors.toList());

    @Autowired
    public ProteinService(ProteinRepository repository) {
        this.repository = repository;
    }

    @PostConstruct
    public void activate() throws IOException {
        allProteins = Stream.of("1m0l").collect(Collectors.toList());

        // clear database
//        repository.deleteAll();
//
//        allProteins.forEach(pdbId -> java.util.concurrent.Executors.newWorkStealingPool().submit(() -> process(pdbId)));
    }

    private void process(String pdbId) {
        try {
            logger.info("fetching information for {}", pdbId);
            Protein protein = OPMDatabaseQuery.parseAnnotatedProteinById(pdbId);
            protein.setTitle(ProteinParser.source(pdbId).parse().getTitle());
            computeFeatures(protein);
            repository.save(new ExplorerProtein(protein));
            logger.info("persisted {}", protein.getName());
        } catch (Exception e) {
            logger.error("gathering information for {} failed: {}", pdbId, e.getLocalizedMessage());
        }
    }

    private void computeFeatures(Protein protein) {
        featureProviders.forEach(featureProvider -> {
            logger.info("processing {} by {}", protein.getName(), featureProvider.getClass().getSimpleName());
            featureProvider.process(protein);
        });
    }

    public List<String> getAllProteins() {
        return allProteins;
    }

    public ExplorerProtein getProtein(String pdbid) {
        List<ExplorerProtein> result = repository.findByTheProteinsName(pdbid);
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
