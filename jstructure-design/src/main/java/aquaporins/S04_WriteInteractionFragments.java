package aquaporins;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.parser.ProteinParser;
import de.bioforscher.jstructure.parser.plip.PLIPInteractionContainer;
import de.bioforscher.jstructure.parser.plip.PLIPParser;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Gather all occurrences of PLIP interactions, align them and write fragments.
 * Created by S on 29.04.2017.
 */
public class S04_WriteInteractionFragments {
    public static void main(String[] args) throws IOException {
        List<PLIPInteractionContainer> containers = Files.list(Paths.get(System.getProperty("user.home") + "/git/phd_sb_repo/data/aquaporins/plip/"))
                .map(S04_WriteInteractionFragments::handlePLIPFile)
                .collect(Collectors.toList());

        containers.forEach(container -> {
            container.getPiStackings();
        });
    }

    private static PLIPInteractionContainer handlePLIPFile(Path path) {
        Chain chain = ProteinParser.source(Paths.get(System.getProperty("user.home") +
                "/git/phd_sb_repo/data/aquaporins/structures/" + path.toFile().getName().split("_")[0] + ".pdb"))
                .parse()
                .select()
                .chainName(path.toFile().getName().split("_")[1].split("\\.")[0])
                .asChain();

        System.out.println(chain.getParentProtein().getName() + "_" + chain.getChainId());

        try {
            return new PLIPInteractionContainer(PLIPParser.parse(chain,
                    Files.lines(path).collect(Collectors.joining(System.lineSeparator()))));
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }
}
