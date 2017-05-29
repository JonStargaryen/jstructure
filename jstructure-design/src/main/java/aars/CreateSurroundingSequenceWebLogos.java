package aars;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.GroupPrototype;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Optional;
import java.util.stream.Collectors;

/**
 * Employ geometry tsvs to generate weblogos.
 * Created by bittrich on 4/19/17.
 */
class CreateSurroundingSequenceWebLogos {
    public static void main(String[] args) throws IOException {
        handleFile("/home/bittrich/git/aars_analysis/data/geometry/bbbrackets_geometry.tsv", 1);
        handleFile("/home/bittrich/git/aars_analysis/data/geometry/argtweezer_geometry.tsv", 2);
    }

    private static void handleFile(String filepath, int classNumber) throws IOException {
        final int[] indices = classNumber == 1 ? new int[] { 326, 1310 } : new int[] { 884, 1765 };

        System.out.println(filepath);

        System.out.println("position 1:");
        Files.lines(Paths.get(filepath))
                .filter(line -> !line.startsWith("id"))
                .map(line -> line.split("\t"))
                .map(split -> {
                    Protein protein = ProteinParser.source(Paths.get("/home/bittrich/git/aars_analysis/data/msa/C" + classNumber + "/renumbered_structures/" + split[0].split("_")[0] + "_renum.pdb")).parse();
                    Chain chain = protein.select().chainName(split[0].split("_")[1]).asChain();
                    int groupIndex = chain.getGroups().indexOf(chain.select().residueNumber(indices[0]).asGroup());
                    return chain.getGroups()
                            .subList(groupIndex - 3, groupIndex + 4)
                            .stream()
                            .map(Group::getGroupPrototype)
                            .map(GroupPrototype::getOneLetterCode)
                            .map(Optional::get)
                            .collect(Collectors.joining());
                })
                .forEach(System.out::println);
        System.out.println();

        System.out.println("position 2:");
        Files.lines(Paths.get(filepath))
                .filter(line -> !line.startsWith("id"))
                .map(line -> line.split("\t"))
                .map(split -> {
                    Protein protein = ProteinParser.source(Paths.get("/home/bittrich/git/aars_analysis/data/msa/C" + classNumber + "/renumbered_structures/" + split[0].split("_")[0] + "_renum.pdb")).parse();
                    Chain chain = protein.select().chainName(split[0].split("_")[1]).asChain();
                    int groupIndex = chain.getGroups().indexOf(chain.select().residueNumber(indices[1]).asGroup());
                    return chain.getGroups()
                            .subList(groupIndex - 3, groupIndex + 4)
                            .stream()
                            .map(Group::getGroupPrototype)
                            .map(GroupPrototype::getOneLetterCode)
                            .map(Optional::get)
                            .collect(Collectors.joining());
                })
                .forEach(System.out::println);
        System.out.println();
    }
}
