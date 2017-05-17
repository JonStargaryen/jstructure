package studies.gmlvq;

import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import de.bioforscher.jstructure.feature.interactions.PLIPAnnotator;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.interactions.PLIPInteraction;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.stream.Collectors;

/**
 * Take the original itemset miner data set and enrich it with PLIP interaction counts.
 * Created by bittrich on 4/18/17.
 */
public class EnrichItemSetMinderDataSetWithInteractionData {
    public static void main(String[] args) throws IOException {
        //TODO global pattern of external files and resources - design should not have resources at all
        String output = Files.lines(Paths.get("/home/bittrich/git/gmlvq_main/data/itemset_miner/PF00127/PF00127.arff"))
                .map(EnrichItemSetMinderDataSetWithInteractionData::handleLine)
                .filter(line -> !line.isEmpty())
                .collect(Collectors.joining(System.lineSeparator()));

        //TODO manually update header
        /*
@ATTRIBUTE halogen1       NUMERIC
@ATTRIBUTE hydrogen1      NUMERIC
@ATTRIBUTE piCation1      NUMERIC
@ATTRIBUTE piStacking1    NUMERIC
@ATTRIBUTE saltBridges1   NUMERIC
@ATTRIBUTE waterBridges1  NUMERIC
         */

        Files.write(Paths.get("/home/bittrich/git/gmlvq_main/data/itemset_miner/PF00127/PF00127_plip.arff"), output.getBytes());
    }

    private static final PLIPAnnotator plipAnnotator = new PLIPAnnotator();

    private static String handleLine(String line) {
        if(line.startsWith("@")) {
            return line;
        }

        try {
            String[] split = line.split(",");
            String[] idSplit = split[0].split("_");

            String pdbId = idSplit[0];

            Protein protein = ProteinParser.source(pdbId).parse();
            plipAnnotator.process(protein);
            PLIPInteractionContainer plipInteractionContainer = protein.getFeatureContainer().getFeature(PLIPInteractionContainer.class);


            Group residue1 = protein.select()
                    .chainName(idSplit[1].split("-")[0])
                    .residueNumber(Integer.valueOf(idSplit[1].split("-")[1].substring(1)))
                    .asGroup();
            String interactionString1 = composeInteractionString(plipInteractionContainer, residue1);
            Group residue2 = protein.select()
                    .chainName(idSplit[2].split("-")[0])
                    .residueNumber(Integer.valueOf(idSplit[2].split("-")[1].substring(1)))
                    .asGroup();
            String interactionString2 = composeInteractionString(plipInteractionContainer, residue2);
            Group residue3 = protein.select()
                    .chainName(idSplit[3].split("-")[0])
                    .residueNumber(Integer.valueOf(idSplit[3].split("-")[1].substring(1)))
                    .asGroup();
            String interactionString3 = composeInteractionString(plipInteractionContainer, residue3);

            String output = split[0] + "," +
                    split[1] + "," +
                    split[2] + "," +
                    split[3] + "," +
                    interactionString1 + "," +
                    split[4] + "," +
                    split[5] + "," +
                    split[6] + "," +
                    interactionString2 + "," +
                    split[7] + "," +
                    split[8] + "," +
                    split[9] + "," +
                    interactionString3 + "," +
                    split[10];
            System.out.println(output);
            return output;
        } catch (Exception e) {
            e.printStackTrace();
            return "";
        }
    }

    private static String composeInteractionString(PLIPInteractionContainer plipInteractionContainer, Group group) {
        long halogenBonds = plipInteractionContainer.getHalogenBonds().stream()
                .filter(halogenBond -> describesGroup(halogenBond, group))
                .count();
        long hydrogenBonds = plipInteractionContainer.getHydrogenBonds().stream()
                .filter(hydrogenBond -> describesGroup(hydrogenBond, group))
                .count();
        long piCation = plipInteractionContainer.getPiCationInteractions().stream()
                .filter(piCationInteraction -> describesGroup(piCationInteraction, group))
                .count();
        long piStacking = plipInteractionContainer.getPiStackings().stream()
                .filter(piStackingInteraction -> describesGroup(piStackingInteraction, group))
                .count();
        long saltBridges = plipInteractionContainer.getSaltBridges().stream()
                .filter(saltBridge -> describesGroup(saltBridge, group))
                .count();
        long waterBridges = plipInteractionContainer.getWaterBridges().stream()
                .filter(waterBridge -> describesGroup(waterBridge, group))
                .count();
        return halogenBonds + "," +
                hydrogenBonds + "," +
                piCation + "," +
                piStacking + "," +
                saltBridges + "," +
                waterBridges;
    }

    private static boolean describesGroup(PLIPInteraction plipInteraction, Group group) {
        return plipInteraction.getPartner1().equals(group) || plipInteraction.getPartner2().equals(group);
    }
}
