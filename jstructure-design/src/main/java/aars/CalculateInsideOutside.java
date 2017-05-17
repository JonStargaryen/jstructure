package aars;

import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceArea;
import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceAreaCalculator;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.family.AminoAcidFamily;
import de.bioforscher.jstructure.parser.ProteinParser;

import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Count inside/outside occurrences of all amino acids in aaRS structures.
 * Created by bittrich on 3/30/17.
 */
class CalculateInsideOutside {
    private static final AccessibleSurfaceAreaCalculator asa = new AccessibleSurfaceAreaCalculator();

    public static void main(String[] args) {
        List<Protein> chains = AARSConstants.lines(Paths.get("/home/bittrich/git/aars_analysis/data/aars_main_table_catalytic.csv"))
                .filter(line -> !line.startsWith("class"))
                .filter(line -> line.startsWith("2"))
                .filter(line -> line.endsWith("1"))
                .map(line -> line.split(","))
                .map(split -> ProteinParser.source(split[2]).parse().select().chainName(split[13]).asChainContainer())
                .peek(System.out::println)
                .map(Protein.class::cast)
                .peek(asa::process)
                .collect(Collectors.toList());

        Map<String, Long> exposedMap = chains.stream()
                .flatMap(Protein::aminoAcids)
                // exposed
                .filter(group -> group.getFeatureContainer().getFeature(AccessibleSurfaceArea.class).getRelativeAccessibleSurfaceArea() > 0.16)
                .collect(Collectors.groupingBy(Group::getThreeLetterCode, Collectors.counting()));

        Map<String, Long> buriedMap = chains.stream()
                .flatMap(Protein::aminoAcids)
                // buried
                .filter(group -> group.getFeatureContainer().getFeature(AccessibleSurfaceArea.class).getRelativeAccessibleSurfaceArea() < 0.16)
                .collect(Collectors.groupingBy(Group::getThreeLetterCode, Collectors.counting()));

        System.out.println("aa\texposed\tburied");
        Stream.of(AminoAcidFamily.values())
                .map(AminoAcidFamily::getThreeLetterCode)
                .map(tlc -> tlc + "\t" + exposedMap.get(tlc) + "\t" + buriedMap.get(tlc))
                .forEach(System.out::println);
    }
}
