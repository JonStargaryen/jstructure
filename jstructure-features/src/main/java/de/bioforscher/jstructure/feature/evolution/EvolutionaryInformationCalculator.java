package de.bioforscher.jstructure.feature.evolution;

import de.bioforscher.jstructure.align.impl.LocalBlastWrapper;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class EvolutionaryInformationCalculator extends FeatureProvider {
    private final LocalBlastWrapper localBlastWrapper;

    public EvolutionaryInformationCalculator() {
        this.localBlastWrapper = new LocalBlastWrapper();
    }

    @Override
    protected void processInternally(Structure protein) {
        protein.chainsWithAminoAcids().forEach(this::processInternally);
    }

    private void processInternally(Chain chain) {
        LocalBlastWrapper.PsiBlastResult psiBlastResult = localBlastWrapper.executePsiBlastUniref50(chain.getAminoAcidSequence());
        assignPsiBlastResultToChain(chain, psiBlastResult);
    }

    public void assignPsiBlastResultToChain(Chain chain, LocalBlastWrapper.PsiBlastResult psiBlastResult) {
        List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());
        for (int i = 0; i < aminoAcids.size(); i++) {
            try {
                aminoAcids.get(i).getFeatureContainer().addFeature(new EvolutionaryInformation(this,
                        psiBlastResult.getExchanges().get(i),
                        psiBlastResult.getInformation().get(i)));
            } catch (IndexOutOfBoundsException e) {
                // happens when no homologs where found in PSI-BLAST run
            }
        }
    }

    public LocalBlastWrapper.PsiBlastResult composePsiBlastResult(Stream<String> outputStream, Stream<String> matrixStream) {
        return localBlastWrapper.composePsiBlastResult(outputStream, matrixStream);
    }
}
