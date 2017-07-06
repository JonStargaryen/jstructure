package studies.poznan;

import org.junit.Test;

/**
 * Test the capabilities of the advanced CLI.
 * Created by bittrich on 2/6/17.
 */
public class AdvancedEnergyProfileRunnerTest {
//    @Test
//    public void computeAndPredictEnergyProfile() {
//        AdvancedEnergyProfileRunner.main(basepath + "1a8h.pdb", basepath + "1a8h_computed.ep");
//        AdvancedEnergyProfileRunner.main(protein1a8h.getAminoAcidSequence(), basepath + "1a8h_predicted.ep");
//        AdvancedEnergyProfileRunner.main(basepath + "1asy.pdb", basepath + "1asy_computed.ep");
//        AdvancedEnergyProfileRunner.main(protein1a8h.getAminoAcidSequence(), basepath + "1asy_predicted.ep");
//    }
//
//    @Test
//    public void alignEnergyProfile() {
//        AdvancedEnergyProfileRunner.main(basepath + "1a8h.pdb", protein1a8h.getAminoAcidSequence(), basepath + "1a8h-egor.epa");
//        AdvancedEnergyProfileRunner.main(basepath + "1asy.pdb", protein1a8h.getAminoAcidSequence(), basepath + "1asy-egor.epa");
//        AdvancedEnergyProfileRunner.main(basepath + "1asy.pdb", basepath + "1asz.pdb", basepath + "1asy-1asz.epa");
//    }

    @Test(expected = IllegalArgumentException.class)
    public void failWithInvalidSequence() {
        AdvancedEnergyProfileRunner.main("KSADQPOWIEPQWOKLASDJKLKSJADLKJQWOPEIOLSKADKLJ", "");
    }

    @Test(expected = IllegalArgumentException.class)
    public void shouldFailDueToNoArguments() {
        AdvancedEnergyProfileRunner.main();
    }

    @Test(expected = IllegalArgumentException.class)
    public void shouldFailDueToTooManyArguments() {
        AdvancedEnergyProfileRunner.main(new String[10]);
    }

    @Test(expected = IllegalArgumentException.class)
    public void shouldFailDueToInvalidFilepath() {
        AdvancedEnergyProfileRunner.main("/invalid/filepath");
    }
}
