package de.bioforscher.jstructure.mutation.impl;

import de.bioforscher.mmm.ItemsetMinerRunner;
import de.bioforscher.mmm.io.DataPointReaderConfiguration;
import de.bioforscher.mmm.model.ItemsetComparatorType;
import de.bioforscher.mmm.model.configurations.ItemsetMinerConfiguration;
import de.bioforscher.mmm.model.configurations.metrics.CohesionMetricConfiguration;
import de.bioforscher.mmm.model.configurations.metrics.ConsensusMetricConfiguration;
import de.bioforscher.mmm.model.configurations.metrics.SupportMetricConfiguration;
import de.bioforscher.mmm.model.enrichment.DataPointEnricherType;
import de.bioforscher.singa.chemistry.physical.model.StructuralEntityFilter;

import java.nio.file.Path;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * Bridge to the molecular miner.
 * Created by bittrich on 7/11/17.
 */
class MolecularMinerBridge {
    private static final ExecutorService executorService = Executors.newWorkStealingPool();

    /**
     * Submit a new job for the molecular miner.
     * @return the instance to-be-created wrapped as Future object
     * @throws ExecutionException thrown upon wrong configuration or internal exceptions
     */
    static Future<ItemsetMinerRunner> submitJob(Path structurePath, Path outputPath) throws ExecutionException {
        try {
            ItemsetMinerConfiguration<String> configuration = createConfiguration(structurePath, outputPath);
            return executorService.submit(() -> new ItemsetMinerRunner(configuration));
        } catch (Exception e) {
            throw new ExecutionException(e);
        }
    }

    private static ItemsetMinerConfiguration<String> createConfiguration(Path structurePath, Path outputPath) {
        ItemsetMinerConfiguration<String> configuration = new ItemsetMinerConfiguration<>();
        configuration.setInputListLocation(null);
        configuration.setInputDirectoryLocation(structurePath.toFile().getAbsolutePath());
        configuration.setOutputLocation(outputPath.toFile().getAbsolutePath());

        DataPointReaderConfiguration dataPointReaderConfiguration = new DataPointReaderConfiguration();
        dataPointReaderConfiguration.setPdbLocation(null);
        dataPointReaderConfiguration.setChainListSeparator("\t");
        dataPointReaderConfiguration.setParseLigands(false);
        dataPointReaderConfiguration.setParseNucleotides(false);
        dataPointReaderConfiguration.setParseWater(false);
        configuration.setDataPointReaderConfiguration(dataPointReaderConfiguration);

        configuration.setDataPointEnricherType(DataPointEnricherType.INTRA_CHAIN_INTERACTION);

        SupportMetricConfiguration<String> simpleMetrics = new SupportMetricConfiguration<>();
        simpleMetrics.setMinimalSupport(0.6);
        configuration.addSimpleMetricConfiguration(simpleMetrics);

        CohesionMetricConfiguration<String> extractionMetric = new CohesionMetricConfiguration<>();
        extractionMetric.setMaximalCohesion(8.0);
        extractionMetric.setVertexOne(false);
        extractionMetric.setLevelOfParallelism(-1);
        configuration.setExtractionMetricConfiguration(extractionMetric);

        ConsensusMetricConfiguration<String> consensus = new ConsensusMetricConfiguration<>();
        consensus.setMaximalConsensus(0.8);
        consensus.setClusterCutoffValue(0.6);
        consensus.setLevelOfParallelism(-1);
        consensus.setAtomFilterType(StructuralEntityFilter.AtomFilterType.ARBITRARY);
        consensus.setRepresentationSchemeType(null);
        consensus.setAlignWithinClusters(true);
        configuration.addExtractionDependentMetricConfiguration(consensus);

        configuration.setItemsetComparatorType(ItemsetComparatorType.CONSENSUS);

        return configuration;
    }
}
