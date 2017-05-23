package studies.gmlvq.cmh;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Element;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.parser.ProteinParser;
import studies.StudyConstants;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.OptionalDouble;
import java.util.stream.Stream;

/**
 * GMLVQ_MAIN misclassifies some instances - check whether functionals are actually non-functional and vice versa.
 * Created by S on 05.05.2017.
 */
class Z01_InvestigateMisclassifiedInstances {
    private static final String FILENAME = StudyConstants.GMLVQ_MAIN + "data/itemset_miner/PF00127/marika/misclassified.txt";

    public static void main(String[] args) throws IOException {
        Files.lines(Paths.get(FILENAME))
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
                            .mapToDouble(coordinates -> atom.algebra().distanceFast(coordinates))
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
