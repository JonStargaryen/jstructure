package plip;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProviderRegistry;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.parser.ProteinParser;
import de.bioforscher.jstructure.parser.plip.PLIPAnnotator;
import de.bioforscher.jstructure.parser.plip.interaction.PLIPInteraction;
import design.DesignConstants;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collection;
import java.util.stream.Collectors;

/**
 * Mine for PLIP interactions in the data set.
 * Created by bittrich on 2/9/17.
 */
@Deprecated
public class S01_MinePLIPInteractions {
    private static final int MINIMAL_SEQUENCE_SEPARATION = 7;
    static final String BASE_PATH = "/home/bittrich/git/phd_sb_repo/data/dataset/nrpdbtm/";
    private static final String DELIMITER = "\t";

    public static void main(String[] args) throws IOException {
        AbstractFeatureProvider plipAnnotator = FeatureProviderRegistry.resolve(PLIPAnnotator.PLIP_INTERACTIONS);

        String listPath = BASE_PATH + "pdbtm_alpha_nr.list.txt";
        String output = Files.lines(Paths.get(listPath))
                .peek(System.out::println)
                .map(line -> line.substring(0, 4))
                .distinct()
                .map(ProteinParser::parseProteinById)
                .map(protein -> {
                    plipAnnotator.process(protein);
                    return protein.getFeatureAsList(PLIPInteraction.class, PLIPAnnotator.PLIP_INTERACTIONS);
                })
                .flatMap(Collection::stream)
                //filter for wide-range interactions
                .filter(plipInteraction -> Math.abs(plipInteraction.getPartner1().getResidueNumber() - plipInteraction.getPartner2().getResidueNumber()) > MINIMAL_SEQUENCE_SEPARATION)
                // filter for residue-residue-interactions only
                .filter(plipInteraction -> plipInteraction.getPartner1().isAminoAcid() && plipInteraction.getPartner2().isAminoAcid())
                .map(plipInteraction -> {
                    Group partner1 = plipInteraction.getPartner1();
                    Group partner2 = plipInteraction.getPartner2();
                    Chain chain = partner1.getParentChain();
                    String pdbId = chain.getParentProtein().getName();
                    String chainId = chain.getChainId();
                    return pdbId + DELIMITER + chainId + DELIMITER + partner1.getThreeLetterCode() + DELIMITER + partner1.getResidueNumber()
                            + DELIMITER + partner2.getThreeLetterCode() + DELIMITER + partner2.getResidueNumber() + DELIMITER + plipInteraction.getClass().getSimpleName();
                })
                .collect(Collectors.joining(System.lineSeparator(),
                        "pdbId" + DELIMITER + "chainId" + DELIMITER + "aa1" + DELIMITER + "resNum1" + DELIMITER + "aa2" + DELIMITER + "resNum2" + DELIMITER + "type" + System.lineSeparator(),
                        ""));

        DesignConstants.write(Paths.get(BASE_PATH + "interactions/plip-data.tsv"), output.getBytes());

                // write to distinct files - redundant data - slow, clogs git
//                .forEach(protein -> {
//                        plipAnnotator.process(protein);
//
//                        List<OldPLIPInteraction> plipInteractions = protein.getFeatureAsList(OldPLIPInteraction.class, PLIPAnnotator.PLIP_INTERACTIONS).stream()
//                                // filter for wide-range interactions
//                                .filter(plipInteraction -> Math.abs(plipInteraction.getPartner1().getResidueNumber() - plipInteraction.getPartner2().getResidueNumber()) > MINIMAL_SEQUENCE_SEPARATION)
//                                // filter for residue-residue-interactions only
//                                .filter(plipInteraction -> plipInteraction.getPartner1().isAminoAcid() && plipInteraction.getPartner2().isAminoAcid())
//                                // filter for non-hydrophobic interactions as they are numerous and potentially only provided little pressure on the arrangement of the protein chain
//            //                                .filter(plipInteraction -> !plipInteraction.getPlipInteractionType().equals(OldPLIPInteractionType.HYDROPHOBIC_INTERACTION))
//                                .collect(Collectors.toList());
//
//                        if(plipInteractions.isEmpty()) {
//                            return;
//                        }
//
//                        String pdbId = plipInteractions.get(0).getPartner1().getParentGroup().getParentGroup().getName();
//                        Path proteinPath = Paths.get(PROJECT_PATH + "interactions/" + pdbId.toLowerCase() + "/");
//                        DesignConstants.makeDirectoryIfAbsent(proteinPath);
//
//                        System.out.println(plipInteractions.stream().collect(Collectors.groupingBy(OldPLIPInteraction::getPlipInteractionType, Collectors.counting())));
//
//                        plipInteractions.forEach(plipInteraction -> {
//                            Group partner1 = plipInteraction.getPartner1();
//                            Group partner2 = plipInteraction.getPartner2();
//                            String chainId = partner1.getParentGroup().getChainId();
//                            DesignConstants.write(proteinPath.resolve(plipInteraction.getPlipInteractionShortName() + "-" + chainId + "-" + partner1.getThreeLetterCode() + partner1.getResidueNumber() + "-" + partner2.getThreeLetterCode() + partner2.getResidueNumber() + ".pdb"),
//                                    plipInteraction.getPartner1().composePDBRecord().concat(plipInteraction.getPartner2().composePDBRecord()).getBytes());
//                        });
//                    });
    }
}
