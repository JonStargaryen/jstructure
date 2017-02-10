package plip;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProviderRegistry;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import de.bioforscher.jstructure.parser.plip.PLIPAnnotator;
import de.bioforscher.jstructure.parser.plip.PLIPInteraction;
import design.DesignConstants;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Mine for PLIP interactions in the data set.
 * Created by bittrich on 2/9/17.
 */
public class S01_MinePLIPInteractions {
    private static final int MINIMAL_SEQUENCE_SEPARATION = 7;
    static final String BASE_PATH = "/home/bittrich/git/phd_sb_repo/data/dataset/nrpdbtm/";

    public static void main(String[] args) throws IOException {
        AbstractFeatureProvider plipAnnotator = FeatureProviderRegistry.resolve(PLIPAnnotator.PLIP_INTERACTIONS);

        String listPath = BASE_PATH + "pdbtm_alpha_nr.list.txt";
        List<Protein> proteins = Files.lines(Paths.get(listPath))
                .peek(System.out::println)
                .map(line -> line.substring(0, 4))
                .distinct()
                .map(ProteinParser::parseProteinById)
                .collect(Collectors.toList());

        proteins.forEach(protein -> {
            plipAnnotator.process(protein);

            List<PLIPInteraction> plipInteractions = protein.getFeatureAsList(PLIPInteraction.class, PLIPAnnotator.PLIP_INTERACTIONS).stream()
                    // filter for wide-range interactions
                    .filter(plipInteraction -> Math.abs(plipInteraction.getPartner1().getResidueNumber() - plipInteraction.getPartner2().getResidueNumber()) > MINIMAL_SEQUENCE_SEPARATION)
                    // filter for residue-residue-interactions only
                    .filter(plipInteraction -> plipInteraction.getPartner1().isAminoAcid() && plipInteraction.getPartner2().isAminoAcid())
                    // filter for non-hydrophobic interactions as they are numerous and potentially only provided little pressure on the arrangement of the protein chain
//                                .filter(plipInteraction -> !plipInteraction.getPlipInteractionType().equals(PLIPInteractionType.HYDROPHOBIC))
                    .collect(Collectors.toList());

            if(plipInteractions.isEmpty()) {
                return;
            }

            String pdbId = plipInteractions.get(0).getPartner1().getParentChain().getParentProtein().getName();
            Path proteinPath = Paths.get(BASE_PATH + "interactions/" + pdbId.toLowerCase() + "/");
            DesignConstants.makeDirectoryIfAbsent(proteinPath);

            System.out.println(plipInteractions.stream().collect(Collectors.groupingBy(PLIPInteraction::getPlipInteractionType, Collectors.counting())));

            plipInteractions.forEach(plipInteraction -> {
                Group partner1 = plipInteraction.getPartner1();
                Group partner2 = plipInteraction.getPartner2();
                String chainId = partner1.getParentChain().getChainId();
                DesignConstants.write(proteinPath.resolve(plipInteraction.getPlipInteractionShortName() + "-" + chainId + "-" + partner1.getThreeLetterCode() + partner1.getResidueNumber() + "-" + partner2.getThreeLetterCode() + partner2.getResidueNumber() + ".pdb"),
                        plipInteraction.getPartner1().composePDBRecord().concat(plipInteraction.getPartner2().composePDBRecord()).getBytes());
            });
        });
    }
}
