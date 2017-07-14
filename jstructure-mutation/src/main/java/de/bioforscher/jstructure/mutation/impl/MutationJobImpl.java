package de.bioforscher.jstructure.mutation.impl;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.mutation.MutationJob;
import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry;

import java.time.LocalDateTime;
import java.util.List;
import java.util.UUID;

/**
 * The implementation of a mutation job.
 * Created by bittrich on 7/14/17.
 */
public class MutationJobImpl implements MutationJob {
    private final String jobName;
    private final UUID uuid;
    private final LocalDateTime creationTime;
    private final Chain originalChain;
    private final Chain referenceChain;
    private List<UniProtEntry> homologousSequences;
    private List<Chain> homologousPdbChains;

    public MutationJobImpl(String jobName,
                           Chain originalChain) {
        this.jobName = jobName;
        this.uuid = UUID.randomUUID();
        this.creationTime = LocalDateTime.now();
        this.originalChain = originalChain;
        // deep-clone original chain - ignore feature map entries - ignore all other chains
        Chain clonedChain = (Chain) originalChain.createCopy();
        Protein protein = new Protein(clonedChain.getChainIdentifier().getProteinIdentifier());
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
    public LocalDateTime getCreationTime() {
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

    public void setHomologousSequences(List<UniProtEntry> homologousSequences) {
        this.homologousSequences = homologousSequences;
    }

    @Override
    public List<UniProtEntry> getHomologousSequences() {
        return homologousSequences;
    }

    public void setHomologousPdbChains(List<Chain> homologousPdbChains) {
        this.homologousPdbChains = homologousPdbChains;
    }

    @Override
    public List<Chain> getHomologousPdbChains() {
        return homologousPdbChains;
    }
}
