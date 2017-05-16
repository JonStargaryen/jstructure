package aquaporins;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Element;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.parser.ProteinParser;
import de.bioforscher.jstructure.parser.plip.PLIPInteractionContainer;
import de.bioforscher.jstructure.parser.plip.PLIPParser;
import de.bioforscher.jstructure.parser.plip.interaction.PLIPInteraction;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Gathers PLIP interaction data and writes the as psuedo-atoms the structures. Aim: visualize aligned structures and
 * screen for conversed respectively shared contacts using PyMOL.
 * Created by bittrich on 4/28/17.
 */
public class S03_WriteStructuresAnnotatedWithPLIPData {
    public static void main(String[] args) throws IOException {
        Files.list(Paths.get(System.getProperty("user.home") + "/git/phd_sb_repo/data/aquaporins/plip/"))
                .forEach(S03_WriteStructuresAnnotatedWithPLIPData::handlePLIPFile);
    }

    private static void handlePLIPFile(Path path) {
        Chain chain = ProteinParser.source(Paths.get(System.getProperty("user.home") +
                "/git/phd_sb_repo/data/aquaporins/structures/" + path.toFile().getName().split("_")[0] + ".pdb"))
                .parse()
                .select()
                .chainName(path.toFile().getName().split("_")[1].split("\\.")[0])
                .asChain();

        System.out.println(chain.getParentProtein().getName() + "_" + chain.getChainId());

        try {
            PLIPInteractionContainer interactionContainer = new PLIPInteractionContainer(PLIPParser.parse(chain, Files.lines(path).collect(Collectors.joining(System.lineSeparator()))));
            addPseudoGroups(chain, interactionContainer);

            Files.write(Paths.get(System.getProperty("user.home") + "/git/phd_sb_repo/data/aquaporins/interaction-structures/" + path.toFile().getName().split("\\.")[0] + ".pdb"),
                    chain.composePDBRecord().getBytes());
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void addPseudoGroups(Chain chain, PLIPInteractionContainer interactionContainer) {
        addPseudoGroups(chain, 9991, "IHA", interactionContainer.getHalogenBonds());
        addPseudoGroups(chain, 9992, "IHB", interactionContainer.getHydrogenBonds());
        addPseudoGroups(chain, 9993, "IMC", interactionContainer.getMetalComplexes());
        addPseudoGroups(chain, 9994, "IPC", interactionContainer.getPiCationInteractions());
        addPseudoGroups(chain, 9995, "IPS", interactionContainer.getPiStackings());
        addPseudoGroups(chain, 9996, "ISB", interactionContainer.getSaltBridges());
        addPseudoGroups(chain, 9997, "IWB", interactionContainer.getWaterBridges());
    }

    private static void addPseudoGroups(Chain chain, int resn, String interactionName, List<? extends PLIPInteraction> interactions) {
        Group group = new Group("UNK", resn);
        interactions.stream()
                // filter only for sequentially separated interactions
                .filter(interaction -> Math.abs(interaction.getPartner1().getResidueNumber() - interaction.getPartner2().getResidueNumber()) > 10)
                .forEach(interaction -> group.addAtom(new Atom(interactionName, 9999, Element.X, interaction.getRepresentation())));
        chain.addGroup(group);
    }
}
