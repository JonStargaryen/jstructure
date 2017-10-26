package de.bioforscher.jstructure.membrane.modularity.pdbtm;

import de.bioforscher.jstructure.align.impl.LocalBlastWrapper;
import de.bioforscher.jstructure.membrane.MembraneConstants;
import de.bioforscher.jstructure.model.structure.StructureParser;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.stream.Collectors;

public class A02_ComputeEvolutionaryInformationProfile {
    private static final Logger logger = LoggerFactory.getLogger(A02_ComputeEvolutionaryInformationProfile.class);

    public static void main(String[] args) {
        LocalBlastWrapper localBlastWrapper = new LocalBlastWrapper();
        MembraneConstants.lines(MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY.resolve("ids.list"))
                .forEach(id -> {
                    logger.info("processing {}", id);
                    String pdbId = id.split("_")[0];
                    String chainId = id.split("_")[1];
                    String sequence = StructureParser.source(pdbId)
                            .minimalParsing(true)
                            .parse()
                            .select()
                            .chainName(chainId)
                            .asChain()
                            .getAminoAcidSequence();
                    logger.info("sequence: {}", sequence);
                    LocalBlastWrapper.PsiBlastResult psiBlastResult = localBlastWrapper.executePsiBlastUniref50(sequence);
                    logger.info("{} aligned sequences", psiBlastResult.getAccessions().size());
                    MembraneConstants.write(MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY
                            .resolve("evol")
                            .resolve(id + ".evol"),
                            psiBlastResult.getInformation()
                            .stream()
                            .map(String::valueOf)
                            .collect(Collectors.joining(System.lineSeparator())));
                });
    }
}
