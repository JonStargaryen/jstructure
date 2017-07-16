package de.bioforscher.jstructure.mutation;

import de.bioforscher.jstructure.model.structure.Chain;
import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry;

import java.time.LocalDate;
import java.util.List;
import java.util.UUID;

/**
 * The data structure holding all information on a mutation job at chain-level. Can be used to predict the effect at
 * individual amino acids, yet all amino acids share this instance so global features do not have to be computed all
 * over again all the time.
 * Created by bittrich on 7/14/17.
 */
public interface MutationJob {
    /**
     * A unique name for this job.
     * @return the name of this job
     */
    String getJobName();

    /**
     *
     * @return
     */
    UUID getUuid();

    /**
     *
     * @return
     */
    LocalDate getCreationTime();

    /**
     * The original chain which will not be manipulated.
     * @return the original chain
     */
    Chain getOriginalChain();

    /**
     * The working copy of the provided chain which will be manipulated.
     * @return
     */
    Chain getReferenceChain();

    /**
     *
     * @return
     */
    List<UniProtEntry> getHomologousSequences();

    /**
     *
     * @return
     */
    List<Chain> getHomologousPdbChains();
}
