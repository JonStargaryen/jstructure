package de.bioforscher.jstructure.mutation;

import de.bioforscher.jstructure.feature.uniprot.homologous.UniProtHomologousEntryContainer;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.util.List;

/**
 * The data structure of predicted mutation effects.
 * Created by bittrich on 7/11/17.
 */
public interface MutationEffectPrediction {
    String getIdentifier();

    String getQuerySequence();

    Protein getQueryProtein();

    Chain getQueryChain();

    void setHomologousEntryContainer(UniProtHomologousEntryContainer feature);

    UniProtHomologousEntryContainer getHomologousEntryContainer();

    void setHomologousPdbChains(List<Chain> homologousPdbChains);

    List<Chain> getHomologousPdbChains();

    boolean predictMutationEffect(int position, AminoAcid.Family targetAminoAcid);
}
