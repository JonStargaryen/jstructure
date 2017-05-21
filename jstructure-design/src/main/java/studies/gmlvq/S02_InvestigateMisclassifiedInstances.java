package studies.gmlvq;

import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Element;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;

import java.util.stream.Stream;

/**
 * GMLVQ misclassifies some instances - check whether functionals are actually non-functional and vice versa.
 * Created by S on 05.05.2017.
 */
public class S02_InvestigateMisclassifiedInstances {
    private static final String MISCLASSIFIED = "24, 1CUR_A-C138_A-H143_A-M148, functional \n"+
            "37, 1ILS_A-C112_A-H117_A-M121, functional \n"+
            "40,   1J5C_A-C83_A-H86_A-M91, functional \n"+
            "56,   1M9W_A-C83_A-H86_A-M91, functional \n"+
            "68,   1PLA_A-C82_A-H85_A-M90, functional \n"+
            "69,   1PLB_A-C82_A-H85_A-M90, functional \n"+
            "70,   1PLC_A-C84_A-H87_A-M92, functional \n"+
            "72,   1PNC_A-C84_A-H87_A-M92, functional \n"+
            "90,   1TU2_A-C89_A-H92_A-M97, functional \n"+
            "106,   2CJ3_A-C89_A-H92_A-M97, functional \n"+
            "129, 2IWE_A-C112_A-H46_A-M121, functional \n"+
            "130,   2JKW_A-C78_A-H81_A-M86, functional \n"+
            "135,   2P80_D-C78_D-H81_D-M86, functional \n"+
            "147,   2W88_A-C89_A-H92_A-M97, functional \n"+
            "184,   4DP4_X-C84_X-H87_X-M92, functional \n"+
            "199, 4KO6_A-C112_A-H117_A-M121, functional \n"+
            "215,   7PCY_A-C84_A-H87_A-M92, functional \n"+
            "217,   9PCY_A-C84_A-H87_A-M92, functional \n"+
            "226,  2EH0_A-C58_A-H57_A-M103, non-functional \n"+
            "245, 5I26_A-C112_A-H117_A-M121, non-functional \n"+
            "257,   1F56_A-C74_A-H79_A-M84, non-functional \n"+
            "260,   4PAZ_A-C78_A-H81_A-M86, non-functional \n"+
            "263, 1OE3_A-C130_A-H139_A-M144, non-functional \n"+
            "286,   2CBP_A-C79_A-H84_A-M89, non-functional \n"+
            "287, 1FWX_A-C565_A-H569_A-M572, non-functional \n"+
            "295, 3X1E_A-C135_A-H143_A-M148, non-functional \n"+
            "310, 5D4H_B-C136_B-H145_B-M150, non-functional \n"+
            "312, 4KNU_A-C123_A-H131_A-M136, non-functional \n"+
            "331, 1OFL_A-C230_A-H198_A-M196, non-functional \n"+
            "373, 1FWX_A-C561_A-H569_A-M572, non-functional \n"+
            "455, 2OC3_A-C131_A-H228_A-M129, non-functional \n"+
            "493, 2YEV_B-C197_B-H205_B-M208, non-functional \n"+
            "517, 3EH5_B-C149_B-H157_B-M160, non-functional \n"+
            "545, 3AG2_B-C196_B-H204_B-M207, non-functional \n"+
            "559, 2ZWA_B-C442_B-H443_B-M441, non-functional \n"+
            "574, 3H2X_A-C227_A-H226_A-M127, non-functional \n"+
            "598, 2VYC_A-C327_A-H363_A-M530, non-functional \n"+
            "726, 1QLE_B-C216_B-H224_B-M227, non-functional \n"+
            "739, 5CZD_A-C266_A-H216_A-M215, non-functional \n"+
            "771, 3H2X_A-C129_A-H226_A-M127, non-functional \n"+
            "888, 3I26_D-C273_D-H195_D-M274, non-functional \n"+
            "977, 3OR2_A-C309_A-H308_A-M313, non-functional \n"+
            "986, 3G5W_A-C147_A-H146_A-M160, non-functional \n"+
            "1016, 1IBJ_C-C417_C-H421_C-M419, non-functional \n"+
            "1049,   1WEU_A-C39_A-H63_A-M49, non-functional \n"+
            "1061, 3X1E_A-C135_A-H134_A-M148, non-functional \n"+
            "1079, 4KNU_A-C123_A-H122_A-M136, non-functional \n"+
            "1086, 5CUS_C-C137_C-H138_C-M139, non-functional \n"+
            "1106, 5D4H_B-C136_B-H135_B-M150, non-functional \n"+
            "1133, 1OE3_A-C130_A-H129_A-M144, non-functional \n"+
            "1156, 3I2T_A-C129_A-H130_A-M131, non-functional \n";

    public static void main(String[] args) {
        Pattern.compile("\n").splitAsStream(MISCLASSIFIED)
                .map(Instance::new)
                .forEach(System.out::println);
    }

    static class Instance {
        Protein protein;
        Group residue1, residue2, residue3;
        OptionalDouble minimalCopperDistance;
        boolean functional, containsCooper;

        Instance(String line) {
            String[] split = line.split(",\\s+")[1].split("_");
            protein = ProteinParser.source(split[0]).parse();

            residue1 = protein.select()
                    .chainName(split[1].split("-")[0])
                    .residueNumber(Integer.valueOf(split[1].split("-")[1].substring(1)))
                    .asGroup();
            residue2 = protein.select()
                    .chainName(split[2].split("-")[0])
                    .residueNumber(Integer.valueOf(split[2].split("-")[1].substring(1)))
                    .asGroup();
            residue3 = protein.select()
                    .chainName(split[3].split("-")[0])
                    .residueNumber(Integer.valueOf(split[3].split("-")[1].substring(1)))
                    .asGroup();

            functional = !line.contains("non-");
            AtomContainer coppers = protein.select().element(Element.Cu).asAtomContainer();
            minimalCopperDistance = coppers.atoms()
                    .mapToDouble(atom -> Stream.of(residue1, residue2, residue3)
                            .map(residue -> residue.select().alphaCarbonAtoms().asAtom())
                            .map(Atom::getCoordinates)
                            .mapToDouble(coordinates -> LinearAlgebra3D.distanceFast(atom.getCoordinates(), coordinates))
                            .average()
                            .getAsDouble())
                    .map(Math::sqrt)
                    .min();

            containsCooper = minimalCopperDistance.isPresent() && minimalCopperDistance.getAsDouble() < 8.0;
        }

        @Override
        public String toString() {
            return "protein=" + protein.getIdentifier() +
                    ", functional=" + functional +
                    ", copperDistance=" + (minimalCopperDistance.isPresent() ? minimalCopperDistance.getAsDouble() : -1) +
                    ", containsCopper=" + containsCooper;
        }
    }
}
