package de.bioforscher.jstructure.feature.uniprot.homologous;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry;
import uk.ac.ebi.kraken.interfaces.uniprot.features.Feature;
import uk.ac.ebi.kraken.interfaces.uniprot.features.FeatureLocation;
import uk.ac.ebi.kraken.interfaces.uniprot.features.MutagenFeature;
import uk.ac.ebi.uniprot.dataservice.client.Client;
import uk.ac.ebi.uniprot.dataservice.client.ServiceFactory;
import uk.ac.ebi.uniprot.dataservice.client.alignment.blast.*;
import uk.ac.ebi.uniprot.dataservice.client.alignment.blast.input.AlignmentCutoffOption;
import uk.ac.ebi.uniprot.dataservice.client.alignment.blast.input.DatabaseOption;
import uk.ac.ebi.uniprot.dataservice.client.alignment.blast.input.ExpectationOption;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Optional;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Annotates a chain with all information available via UniProt entries of (potentially) homologous proteins.
 * Created by bittrich on 7/10/17.
 */
@FeatureProvider(provides = { UniProtHomologousEntryContainer.class, UniProtFeatureContainer.class })
public class UniProtHomologyAnnotator extends AbstractFeatureProvider {
    private static final Logger logger  = LoggerFactory.getLogger(UniProtHomologyAnnotator.class);
    final UniProtBlastService uniProtBlastService;

    public UniProtHomologyAnnotator() {
        ServiceFactory serviceFactoryInstance = Client.getServiceFactoryInstance();
        uniProtBlastService = serviceFactoryInstance.getUniProtBlastService();
        uniProtBlastService.start();
    }

    @Override
    protected void processInternally(Protein protein) {
        //TODO statelessness and co
        //TODO start/stop of service
        protein.chainsWithAminoAcids()
                .forEach(this::processInternally);
    }

    private void processInternally(Chain chain) {
        String sequence = chain.getAminoAcidSequence();
        List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());
        // initialize all residues with empty container
        aminoAcids.stream()
                .map(AminoAcid::getFeatureContainer)
                .forEach(featureContainer -> featureContainer.addFeature(new UniProtFeatureContainer(this)));

        logger.info("submitting blast query");
        logger.debug("sequence{}>{}{}{}", System.lineSeparator(), chain.getChainIdentifier(), System.lineSeparator(), sequence);
        List<UniProtHit> uniProtHits = runUniProtBlastService(sequence);
        UniProtHomologousEntryContainer uniProtHomologousEntryContainer = new UniProtHomologousEntryContainer(this,
                uniProtHits);
        // set mapping on chain level
        chain.getFeatureContainer().addFeature(uniProtHomologousEntryContainer);

        for (UniProtHit uniProtHit : uniProtHits) {
            UniProtEntry uniProtEntry = uniProtHit.getEntry();
            String accession = uniProtEntry.getPrimaryUniProtAccession().getValue();
            logger.info("processing hit {}", accession);
            Collection<Feature> features = uniProtEntry.getFeatures();

            Alignment alignment = uniProtHit.getSummary().getAlignments().get(0);
            int querySequenceOffset = alignment.getStartQuerySeq() - 1;
            int matchSequenceOffset = alignment.getStartMatchSeq() - 1;
            for (Feature feature : features) {
                if (!(feature instanceof MutagenFeature)) {
                    continue;
                }

                FeatureLocation featureLocation = feature.getFeatureLocation();
                // nothing to do if start is unknown
                if (!featureLocation.isStartAvailable()) {
                    continue;
                }

                int start = featureLocation.getStart() - 1;
                int end = (featureLocation.isEndAvailable() ? featureLocation.getEnd() : start);
                try {
                    IntStream.range(start, end)
                            // respect unaligned, terminal regions of match sequence
                            .map(position -> position - matchSequenceOffset)
                            // map from match to query
                            .mapToObj(position -> mapFromMatchToQuery(alignment, position))
                            .filter(Optional::isPresent)
                            .mapToInt(Optional::get)
                            // respect unaligned, terminal regions of query sequence
                            .map(position -> position + querySequenceOffset)
                            .mapToObj(position -> getSafely(aminoAcids, position))
                            .filter(Optional::isPresent)
                            .map(Optional::get)
                            .map(AminoAcid::getFeatureContainer)
                            .forEach(featureContainer -> featureContainer.getFeature(UniProtFeatureContainer.class)
                                    .addFeature(accession, feature));
                } catch (Exception e) {
                    // happens on regions not covered or present in the SEQRES but does not actually provide coordinates
                }
            }
        }
    }

    private Optional<AminoAcid> getSafely(List<AminoAcid> aminoAcids, int index) {
        try {
            return Optional.of(aminoAcids.get(index));
        } catch (ArrayIndexOutOfBoundsException e) {
            return Optional.empty();
        }
    }

    private Optional<Integer> mapFromMatchToQuery(Alignment alignment, int position) {
        try {
            // operate on partial string, ignore last position for now
            String querySequence = alignment.getQuerySeq().substring(0, position);
            int queryGaps = countGaps(querySequence);
            String matchSequence = alignment.getMatchSeq().substring(0, position);
            int matchGaps = countGaps(matchSequence);

            char charInQuery = alignment.getQuerySeq().charAt(position);
            if (charInQuery == '-') {
                return Optional.empty();
            } else {
                return Optional.of(position + matchGaps - queryGaps);
            }
        } catch (StringIndexOutOfBoundsException e) {
            // happens on regions not covered or present in the SEQRES but does not actually provide coordinates
            return Optional.empty();
        }
    }

    private int countGaps(String sequence) {
        int count = 0;
        for(int i = 0; i < sequence.length(); i++) {
            if(sequence.charAt(i) == '-') {
                count++;
            }
        }
        return count;
    }

    List<UniProtHit> runUniProtBlastService(String querySequence) {
        return runUniProtBlastService(querySequence,
                DatabaseOption.SWISSPROT,
                //TODO ponder on cutoff
                AlignmentCutoffOption.ONE_THOUSAND,
                ExpectationOption.ONE_EXPONENT_MINUS_FOUR);
    }

    private List<UniProtHit> runUniProtBlastService(String querySequence,
                                                   DatabaseOption databaseOption,
                                                   AlignmentCutoffOption alignmentCutoffOption,
                                                   ExpectationOption expectationOption) {
        List<UniProtHit> hits = new ArrayList<>();

        BlastInput input = new BlastInput.Builder(databaseOption, querySequence)
                .withMaximumNumberOfAlignments(alignmentCutoffOption)
                .withExpectation(expectationOption)
                .build();

        CompletableFuture<BlastResult<UniProtHit>> resultFuture = uniProtBlastService.runBlast(input);

        try {
            // blocks until result arrives
            BlastResult<UniProtHit> blastResult = resultFuture.get();
            logger.info("Number of blast hits: " + blastResult.getNumberOfHits());
            // populate result list
            blastResult.hits().forEach(hits::add);
        } catch (ExecutionException e) {
            logger.error(e.getCause().getMessage());
        } catch (InterruptedException e) {
            logger.error(e.getMessage());
        }

        return hits;
    }
}
