package de.bioforscher.explorer.membrane;

import de.bioforscher.explorer.membrane.model.ExplorerChain;
import de.bioforscher.explorer.membrane.model.ModelFactory;
import de.bioforscher.explorer.membrane.model.MultiSequenceAlignment;
import de.bioforscher.explorer.membrane.repository.AlignmentRepository;
import de.bioforscher.explorer.membrane.repository.ChainRepository;
import de.bioforscher.jstructure.parser.pdb.PDBDatabaseQuery;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import javax.annotation.PostConstruct;
import java.io.IOException;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Access to the repositories.
 * Created by bittrich on 3/17/17.
 */
@Service
public class ChainService {
    private static final Logger logger = LoggerFactory.getLogger(ChainService.class);

    private ChainRepository chainRepository;
    private AlignmentRepository alignmentRepository;

    private List<String> allRepresentativeChainIds;

    @Autowired
    public ChainService(ChainRepository chainRepository, AlignmentRepository alignmentRepository) {
        this.chainRepository = chainRepository;
        this.alignmentRepository = alignmentRepository;
    }

    @PostConstruct
    public void activate() throws IOException {
        logger.info("starting protein service");
        allRepresentativeChainIds = Stream.of("1m0l/A").collect(Collectors.toList());

        // clear database
        logger.info("dropping database");
        chainRepository.deleteAll();

        logger.info("initializing database of {} representative protein chains", allRepresentativeChainIds.size());
        allRepresentativeChainIds.forEach(pdbId ->
//                java.util.concurrent.Executors.newWorkStealingPool().submit(() ->
                        processSequenceCluster(pdbId)
//                )
        );
    }

    private void processSequenceCluster(String clusterRepresentativeChainId) {
        logger.info("[{}] spawned thread on representative chain id", clusterRepresentativeChainId);
        String split[] = clusterRepresentativeChainId.split("/");
        List<String> clusterChainIds = PDBDatabaseQuery.fetchSequenceCluster(split[0], split[1]);
        logger.info("[{}] mapped homologous {}", clusterRepresentativeChainId, clusterChainIds);

        List<ExplorerChain> explorerChains = ModelFactory.createChains(clusterRepresentativeChainId, clusterChainIds);

        logger.info("[{}] creating multi-sequence alignment by clustal omega", clusterRepresentativeChainId);
        //TODO case of no homologous
        MultiSequenceAlignment multiSequenceAlignment = ModelFactory.createMultiSequenceAlignment(clusterRepresentativeChainId, explorerChains);

        logger.info("[{}] persisting chains and alignment", clusterRepresentativeChainId);
        chainRepository.save(explorerChains);
        alignmentRepository.save(multiSequenceAlignment);
    }

    public List<String> getAllRepresentativeChainIds() {
        return allRepresentativeChainIds;
    }

    public ExplorerChain getChain(String rawId) {
        System.out.println(rawId);
        // ensure pdbId is always in upper case
        String id = rawId.split("\\.")[0].toUpperCase() + "." + rawId.split("\\.")[1];
        List<ExplorerChain> result = chainRepository.findChainById(id);
        if(result.size() > 0) {
            return result.get(0);
        }

        throw new NoSuchElementException("could not retrieve entry for " + rawId);
    }

    public List<String> getAllChainIds() {
        return chainRepository.findAll().stream()
                .map(ExplorerChain::getId)
                .collect(Collectors.toList());
    }
}
