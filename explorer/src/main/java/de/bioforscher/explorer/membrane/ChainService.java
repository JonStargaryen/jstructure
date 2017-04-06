package de.bioforscher.explorer.membrane;

import de.bioforscher.explorer.membrane.repository.AlignmentRepository;
import de.bioforscher.explorer.membrane.repository.ChainRepository;
import de.bioforscher.explorer.model.ExplorerAlignment;
import de.bioforscher.explorer.model.ExplorerChain;
import de.bioforscher.explorer.model.ExplorerCluster;
import de.bioforscher.explorer.model.ModelFactory;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import javax.annotation.PostConstruct;
import java.io.IOException;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.concurrent.Executors;
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
        allRepresentativeChainIds = Stream.of("1ea5_A").collect(Collectors.toList());

        reinitialize(false);
    }

    private void reinitialize(boolean parallel) {
        // clear database
        logger.info("dropping database");
        chainRepository.deleteAll();
        alignmentRepository.deleteAll();

        logger.info("initializing database of {} representative protein chains", allRepresentativeChainIds.size());

        if(parallel) {
            allRepresentativeChainIds.forEach(pdbId -> Executors.newWorkStealingPool().submit(() -> processSequenceCluster(pdbId)));
        } else {
            allRepresentativeChainIds.forEach(this::processSequenceCluster);
        }
    }

    private void processSequenceCluster(String clusterEntryId) {
        logger.info("[{}] spawned thread on representative chain id", clusterEntryId);
        ExplorerCluster explorerCluster = ModelFactory.createCluster(clusterEntryId);

        logger.info("[{}] persisting chains and alignment", clusterEntryId);
        chainRepository.save(explorerCluster.getExplorerChains());
        alignmentRepository.save(explorerCluster.getExplorerAlignment());
    }

    public List<String> getAllRepresentativeChainIds() {
        return allRepresentativeChainIds;
    }

    public ExplorerChain getChain(String rawId) {
        String[] split = rawId.split("_");
        // ensure pdbId is always in upper case
        String id = split[0].toLowerCase() + "_" + split[1];
        List<ExplorerChain> result = chainRepository.findChainById(id);
        if(result.size() > 0) {
            return result.get(0);
        }

        throw new NoSuchElementException("could not retrieve entry for " + rawId);
    }

    public ExplorerAlignment getAlignment(String representativeId) {
        List<ExplorerAlignment> result = alignmentRepository.findAlignmentByRepresentativeId(representativeId);
        if(result.size() > 0) {
            return result.get(0);
        }

        throw new NoSuchElementException("could not retrieve multi-sequence alignment for representative " + representativeId);
    }

    public List<String> getAllChainIds() {
        return chainRepository.findAll().stream()
                .map(ExplorerChain::getId)
                .collect(Collectors.toList());
    }
}
