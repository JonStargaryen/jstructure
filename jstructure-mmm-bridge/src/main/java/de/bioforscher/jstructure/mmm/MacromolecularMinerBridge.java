package de.bioforscher.jstructure.mmm;

import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.mmm.ItemsetMinerRunner;
import de.bioforscher.mmm.model.configurations.ItemsetMinerConfiguration;

import java.nio.file.Path;
import java.util.concurrent.CompletableFuture;

/**
 * Integration of the macromolecular miner.
 * Created by bittrich on 7/13/17.
 */
public interface MacromolecularMinerBridge {
    /**
     * Submit a new job for the molecular miner.
     * @param structurePath the path were the structure ensemble resides
     * @param outputPath the desired output path
     * @return the instance to-be-created wrapped as Future object
     */
    default CompletableFuture<ItemsetMinerRunner> submitStandardJob(Path structurePath, Path outputPath) {
        return submitJob(getStandardConfiguration(structurePath, outputPath));
    }

    /**
     * @return the standard configuration
     */
    ItemsetMinerConfiguration<String> getStandardConfiguration(Path structurePath, Path outputPath);

    /**
     * @return a dedicated configuration employing Gutteridge-mapping and enforcing a sequence separation value
     */
    ItemsetMinerConfiguration<String> getConservationProfileConfiguration(Path structurePath, Path outputPath);

    /**
     * Submit a new job for the molecular miner.
     * @param itemsetMinerConfiguration the configuration to use
     * @return the instance to-be-created wrapped as Future object
     */
    CompletableFuture<ItemsetMinerRunner> submitJob(ItemsetMinerConfiguration<String> itemsetMinerConfiguration);

    CompletableFuture<Structure> getConservationProfile(Path structurePath, Structure referenceProtein);
}
