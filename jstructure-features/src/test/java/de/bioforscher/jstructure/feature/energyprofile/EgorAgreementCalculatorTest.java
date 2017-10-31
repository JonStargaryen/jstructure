package de.bioforscher.jstructure.feature.energyprofile;

import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.testutil.TestUtils;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

public class EgorAgreementCalculatorTest {
    private Structure structure;
    private AminoAcid aminoAcid;

    @Before
    public void setup() {
        structure = StructureParser.source(TestUtils.getProteinInputStream(TestUtils.SupportedProtein.PDB_1ACJ))
                .minimalParsing(true)
                .parse();
        aminoAcid = structure.aminoAcids()
                .findFirst()
                .get();
    }

    @Test
    public void shouldResolveRequirements() {
        EgorAgreement egorAgreement = aminoAcid.getFeature(EgorAgreement.class);
        Assert.assertNotEquals("computation and prediction match perfectly - probably the correct feature providers where not resolved",
                egorAgreement.getSolvationEnergy(),
                egorAgreement.getEgorPrediction(),
                TestUtils.TOLERANT_ERROR_MARGIN);
    }
}