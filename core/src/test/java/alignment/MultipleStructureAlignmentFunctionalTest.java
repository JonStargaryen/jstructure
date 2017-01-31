package alignment;

import de.bioforscher.jstructure.alignment.multiple.MultipleStructureAlignment;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Run the binding site alignment.
 * Created by bittrich on 1/25/17.
 */
public class MultipleStructureAlignmentFunctionalTest {
    private static final Logger logger = LoggerFactory.getLogger(MultipleStructureAlignmentFunctionalTest.class);

    @Test
    public void shouldPerformMultipleStructureAlignment() throws IOException {
        int arginine1 = 884;
        int arginine2 = 1765;

        logger.info("loading structures...");
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

        MultipleStructureAlignment multipleStructureAlignment = new MultipleStructureAlignment();
        multipleStructureAlignment.align(chains, arginine1, arginine2);

        Map<GroupContainer, GroupContainer> alignedContainers = multipleStructureAlignment.getAlignedContainerMap();

        alignedContainers.entrySet().forEach(entry -> {
            try {
                GroupContainer alignedFullStructure = entry.getKey();
                GroupContainer alignedCore = entry.getValue();
                int numberOfAlignedGroups = alignedCore.getGroups().size();
                Files.write(Paths.get(homePath + "/multiple-alignment/core/" + numberOfAlignedGroups + "-" + alignedCore.getIdentifier() + ".pdb"), alignedCore.composePDBRecord().getBytes());
                Files.write(Paths.get(homePath + "/multiple-alignment/full/" + alignedFullStructure.getIdentifier() + ".pdb"), alignedFullStructure.composePDBRecord().getBytes());
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
        });
    }
}
