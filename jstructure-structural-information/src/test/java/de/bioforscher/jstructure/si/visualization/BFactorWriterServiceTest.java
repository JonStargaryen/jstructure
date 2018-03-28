package de.bioforscher.jstructure.si.visualization;

import de.bioforscher.jstructure.efr.model.Start2FoldResidueAnnotation;
import de.bioforscher.jstructure.efr.parser.Start2FoldXmlParser;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.testutil.FileUtils;
import org.junit.Ignore;
import org.junit.Test;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class BFactorWriterServiceTest {
    @Test
    @Ignore
    public void shouldWriteBFactorAnnotatedFiles() {
        writeFiles("1bdd", "STF0045", Paths.get("/home/bittrich/tmp/"), 185);
        writeFiles("1omu", "STF0016", Paths.get("/home/bittrich/tmp/"), 76);
    }

    private void writeFiles(String pdbId,
                            String entryId,
                            Path outDir,
                            int... experimentIds) {
        Chain chain = StructureParser.fromPdbId(pdbId)
                .parse()
                .getFirstChain();
        List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());

        Start2FoldXmlParser.parseSpecificExperiment(chain,
                FileUtils.newInputStream(Paths.get("/home/bittrich/git/phd_sb_repo/data/start2fold/xml/" + entryId + ".xml")),
                IntStream.of(experimentIds).boxed().collect(Collectors.toList()));

        List<AminoAcid> earlyFoldingResidues = chain.aminoAcids()
                .filter(aminoAcid -> aminoAcid.getFeature(Start2FoldResidueAnnotation.class).isEarly())
                .collect(Collectors.toList());

        List<ContactStructuralInformation> contactStructuralInformation = StructuralInformationParserService.getInstance()
                .parseContactStructuralInformationFile(FileUtils.newInputStream(Paths.get("/home/bittrich/git/phd_sb_repo/data/si/raw/" + entryId + ".out")),
                        earlyFoldingResidues);
        List<ResidueStructuralInformation> residueStructuralInformation = StructuralInformationParserService.getInstance()
                .composeResidueStructuralInformation(aminoAcids,
                        earlyFoldingResidues,
                        contactStructuralInformation);

        BFactorWriterService.getInstance().writePDBFileWithBFactors(chain,
                residueStructuralInformation,
                outDir.resolve(entryId + "_avg.pdb"),
                outDir.resolve(entryId + "_max.pdb"));
    }
}