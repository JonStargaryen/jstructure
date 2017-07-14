package de.bioforscher.jstructure.mmm.impl;

import de.bioforscher.jstructure.mmm.MacromolecularMinerBridge;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.mmm.ItemsetMiner;
import de.bioforscher.mmm.ItemsetMinerRunner;
import de.bioforscher.mmm.io.DataPointReaderConfiguration;
import de.bioforscher.mmm.model.ItemsetComparatorType;
import de.bioforscher.mmm.model.configurations.ItemsetMinerConfiguration;
import de.bioforscher.mmm.model.configurations.metrics.CohesionMetricConfiguration;
import de.bioforscher.mmm.model.configurations.metrics.ConsensusMetricConfiguration;
import de.bioforscher.mmm.model.configurations.metrics.SeparationMetricConfiguration;
import de.bioforscher.mmm.model.configurations.metrics.SupportMetricConfiguration;
import de.bioforscher.mmm.model.enrichment.IntraChainInteractionEnricher;
import de.bioforscher.mmm.model.mapping.rules.ChemicalGroupsMappingRule;
import de.bioforscher.singa.chemistry.physical.model.StructuralEntityFilter;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.concurrent.CompletableFuture;

/**
 * Implementation of the macromolecular miner bridge.
 * Created by bittrich on 7/13/17.
 */
public class MacromolecularMinerBridgeImpl implements MacromolecularMinerBridge {
    private final StructureConversationCalculatorImpl structureConservationCalculator;

    public MacromolecularMinerBridgeImpl() {
        this.structureConservationCalculator = new StructureConversationCalculatorImpl();
    }

    @Override
    public ItemsetMinerConfiguration<String> getStandardConfiguration(Path structurePath, Path outputPath) {
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

        configuration.setDataPointEnricher(new IntraChainInteractionEnricher());

        SupportMetricConfiguration<String> simpleMetrics = new SupportMetricConfiguration<>();
        simpleMetrics.setMinimalSupport(0.8);
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

    @Override
    public ItemsetMinerConfiguration<String> getConservationProfileConfiguration(Path structurePath, Path outputPath) {
        {
            ItemsetMinerConfiguration<String> configuration = getStandardConfiguration(structurePath, outputPath);

            // mine on functional groups
//        configuration.setMappingRule(new FunctionalGroupsMappingRule());
            // mine on gutteridge-mapped amino acids
            configuration.setMappingRule(new ChemicalGroupsMappingRule());

            // force a certain sequence separation
            SeparationMetricConfiguration<String> separationMetricConfiguration = new SeparationMetricConfiguration<>();
            separationMetricConfiguration.setMaximalSeparation(50);
            separationMetricConfiguration.setOptimalSeparation(5);
            separationMetricConfiguration.setMorseWellDepth(500);
            separationMetricConfiguration.setMorseShape(0.2);
            configuration.addExtractionDependentMetricConfiguration(separationMetricConfiguration);

            //TODO a way to change the minimal support

            return configuration;
        }
    }

    @Override
    public CompletableFuture<ItemsetMinerRunner> submitJob(ItemsetMinerConfiguration<String> itemsetMinerConfiguration) {
        CompletableFuture<ItemsetMinerRunner> future = new CompletableFuture<>();
        CompletableFuture.runAsync(() -> {
            try {
                future.complete(new ItemsetMinerRunner(itemsetMinerConfiguration));
            } catch (Exception e) {
                future.completeExceptionally(e);
            }
        });
        return future;
    }

    @Override
    public CompletableFuture<Protein> getConservationProfile(Path structurePath, Protein referenceProtein) {
        try {
            Path outputPath = Files.createTempDirectory("mmm-out");
            return submitJob(getConservationProfileConfiguration(structurePath, outputPath))
                    .thenApply(ItemsetMinerRunner::getItemsetMiner)
                    .thenApply(ItemsetMiner::getTotalExtractedItemsets)
                    .thenApply(extractedItemsets -> {
                        structureConservationCalculator.extractConservationProfile(extractedItemsets, referenceProtein);
                        return referenceProtein;
                    });
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }
}
