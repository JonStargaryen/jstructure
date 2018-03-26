package de.bioforscher.jstructure.si.visualization;

import de.bioforscher.jstructure.efr.model.Start2FoldResidueAnnotation;
import de.bioforscher.jstructure.efr.parser.Start2FoldXmlParser;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.junit.Ignore;
import org.junit.Test;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class StructuralInformationParserServiceTest {
    @Test
    @Ignore
    public void shouldParseStructuralInformationFile() {
        Chain chain = StructureParser.fromPdbId("1bdd")
                .parse()
                .getFirstChain();

        Start2FoldXmlParser.parseSpecificExperiment(chain,
                Paths.get("/home/bittrich/git/phd_sb_repo/data/start2fold/xml/STF0045.xml"),
                Stream.of(185).collect(Collectors.toList()));

        List<AminoAcid> earlyFoldingResidues = chain.aminoAcids()
                .filter(aminoAcid -> aminoAcid.getFeature(Start2FoldResidueAnnotation.class).isEarly())
                .collect(Collectors.toList());

        Path path = Paths.get("/home/bittrich/git/phd_sb_repo/data/si/raw/STF0045.out");
        StructuralInformationParserService.getInstance().parseStructuralInformationFile(path,
                earlyFoldingResidues)
                .forEach(System.out::println);
    }
}