package de.bioforscher.jstructure.mutation.old;

import de.bioforscher.jstructure.feature.uniprot.homologous.UniProtHomologousEntryContainer;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.mutation.old.impl.SequenceConservationProfile;

import java.util.List;
import java.util.Map;

/**
 * The data structure of predicted mutation effects.
 * Created by bittrich on 7/11/17.
 */
public interface MutationJob {
    String getIdentifier();

    String getQuerySequence();

    Structure getQueryProtein();

    Chain getQueryChain();

    void setHomologousEntryContainer(UniProtHomologousEntryContainer feature);

    UniProtHomologousEntryContainer getHomologousEntryContainer();

    void setHomologousPdbChains(List<Chain> homologousPdbChains);

    List<Chain> getHomologousPdbChains();

    Structure getReferenceProtein();

    void setReferenceProtein(Structure referenceProtein);

    Chain getReferenceChain();

    void setReferenceChain(Chain referenceChain);

    Map<String, String> getAlignmentMap();

    void setAlignmentMap(Map<String, String> alignmentMap);

    List<SequenceConservationProfile> getSequenceConservationProfile();

    MutationDescriptor composeMutationDescriptor(int position, AminoAcid.Family targetAminoAcid);

    boolean predictMutationEffect(int position, AminoAcid.Family targetAminoAcid);
}
