package de.bioforscher.jstructure.membrane.modularity.pdbtm;

import de.bioforscher.jstructure.align.impl.LocalBlastWrapper;
import de.bioforscher.jstructure.membrane.MembraneConstants;
import de.bioforscher.jstructure.model.structure.StructureParser;

import java.util.stream.Collectors;

public class B02_ComputeEvolutionaryInformationProfile {
    public static void main(String[] args) {
        LocalBlastWrapper localBlastWrapper = new LocalBlastWrapper();
        MembraneConstants.lines(MembraneConstants.PDBTM_NR_BETA_DATASET_DIRECTORY.resolve("ids.list"))
                .forEach(id -> {
                    System.out.println("processing " + id);
                    String pdbId = id.split("_")[0];
                    String chainId = id.split("_")[1];
                    String sequence = StructureParser.source(pdbId)
                            .minimalParsing(true)
                            .parse()
                            .select()
                            .chainName(chainId)
                            .asChain()
                            .getAminoAcidSequence();
                    System.out.println("sequence: " + sequence);
                    LocalBlastWrapper.PsiBlastResult psiBlastResult = localBlastWrapper.executePsiBlastUniref50(sequence);
                    System.out.println(psiBlastResult.getAccessions().size() + " aligned sequences");
                    MembraneConstants.write(MembraneConstants.PDBTM_NR_BETA_DATASET_DIRECTORY
                            .resolve("evol")
                            .resolve(id + ".evol"),
                            psiBlastResult.getInformation()
                            .stream()
                            .map(String::valueOf)
                            .collect(Collectors.joining(System.lineSeparator())));
                });
    }
}
