package de.bioforscher.jstructure.mutation.impl;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.mutation.MutationJob;
import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry;

import java.time.LocalDate;
import java.util.List;
import java.util.Map;
import java.util.UUID;

/**
 * The implementation of a mutation job.
 * Created by bittrich on 7/14/17.
 */
public class MutationJobImpl implements MutationJob {
    private final String jobName;
    private final UUID uuid;
    private final LocalDate creationTime;
    private final Chain originalChain;
    private final Chain referenceChain;
    private List<UniProtEntry> homologousSequences;
    private List<Chain> homologousPdbChains;
    private Map<String, String> alignmentMap;

    public MutationJobImpl(String jobName,
                           Chain originalChain) {
        this.jobName = jobName;
        this.uuid = UUID.randomUUID();
        this.creationTime = LocalDate.now();
        this.originalChain = originalChain;
        // deep-clone original chain - ignore feature map entries - ignore all other chains
        Chain clonedChain = originalChain.createDeepCopy();
        Structure protein = new Structure(clonedChain.getChainIdentifier().getProteinIdentifier());
        protein.addChain(clonedChain);
        this.referenceChain = clonedChain;
    }

    @Override
    public String getJobName() {
        return jobName;
    }

    @Override
    public UUID getUuid() {
        return uuid;
    }

    @Override
    public LocalDate getCreationTime() {
        return creationTime;
    }

    @Override
    public Chain getOriginalChain() {
        return originalChain;
    }

    @Override
    public Chain getReferenceChain() {
        return referenceChain;
    }

    void setHomologousSequences(List<UniProtEntry> homologousSequences) {
        this.homologousSequences = homologousSequences;
    }

    @Override
    public List<UniProtEntry> getHomologousSequences() {
        return homologousSequences;
    }

    void setHomologousPdbChains(List<Chain> homologousPdbChains) {
        this.homologousPdbChains = homologousPdbChains;
    }

    @Override
    public List<Chain> getHomologousPdbChains() {
        return homologousPdbChains;
    }

    public void setAlignmentMap(Map<String, String> alignmentMap) {
        this.alignmentMap = alignmentMap;
    }

    @Override
    public Map<String, String> getAlignmentMap() {
        return alignmentMap;
    }
}
