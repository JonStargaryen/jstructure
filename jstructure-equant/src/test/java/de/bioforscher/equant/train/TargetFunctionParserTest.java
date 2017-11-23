package de.bioforscher.equant.train;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.testutil.TestUtils;
import org.junit.Test;

public class TargetFunctionParserTest {
    @Test
    public void shouldParseResultFile() {
        Structure structure = StructureParser.source(TestUtils.getResourceAsInputStream("T0540TS001_1")).parse();
        new TargetFunctionParser().process(structure, TestUtils.getResourceAsInputStream("T0540TS001_1.lga"));
        structure.aminoAcids()
                .map(aminoAcid -> aminoAcid.getResidueIdentifier() + " " +
                        StandardFormat.format(aminoAcid.getFeature(TargetFunctionParser.TargetFunction.class).getDistance()) + " " +
                        StandardFormat.format(aminoAcid.getFeature(TargetFunctionParser.TargetFunction.class).getSscore()))
                .forEach(System.out::println);
    }
}