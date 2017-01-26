package alignment;

import de.bioforscher.jstructure.alignment.multiple.MultipleStructureAlignment;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Test;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Run the binding site alignment.
 * Created by bittrich on 1/25/17.
 */
public class BindingSiteAlignmentFunctionalTest {
    @Test
    public void shouldPerformMultipleStructureAlignment() throws IOException {
        int arginine1 = 884;
        int arginine2 = 1765;

        String homePath = System.getProperty("user.home");
        String basePath = homePath + "/git/aars_analysis/data/";
        List<GroupContainer> chains = Files.lines(Paths.get(basePath + "geometry/argtweezer_geometry.tsv"))
                .filter(line -> !line.startsWith("id"))
                .map(line -> line.split("\t"))
                .map(split -> split[0])
                .map(id -> {
                    //TODO selection-API should resolve/infer container name if possible
                    GroupContainer container = Selection.on(ProteinParser.parsePDBFile(basePath + "msa/C2/renumbered_structures/" + id.split("_")[0].toLowerCase() + "_renum.pdb"))
                            .chainName(id.split("_")[1])
                            .asGroupContainer();
                    container.setIdentifier(id);
                    return container;
                })
                .collect(Collectors.toList());

        List<GroupContainer> alignedFragments = new MultipleStructureAlignment().align(chains, arginine1, arginine2);
        for(GroupContainer container : alignedFragments) {
            Files.write(Paths.get(homePath + "/multiple-alignment/" + container.getIdentifier() + ".pdb"), container.composePDBRecord().getBytes());
        }
    }
}
