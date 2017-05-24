package aars;

import de.bioforscher.jstructure.alignment.SVDSuperimposer;
import de.bioforscher.jstructure.alignment.StructureAlignmentResult;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Element;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static aars.AARSConstants.MAIN_TABLE_CATALYTIC_PATH;

/**
 * Aligns binding sites with respect to a present ligand.
 * Created by bittrich on 2/7/17.
 */
@Deprecated
class BindingSiteAligner {
    private static final Logger logger = LoggerFactory.getLogger(BindingSiteAligner.class);

    public static void main(String[] args) throws IOException {
        // load proteins
        Map<String, Protein> proteinMap = Files.lines(Paths.get(MAIN_TABLE_CATALYTIC_PATH))
                // skip header line
                .filter(line -> !line.startsWith("class"))
                // filter for rep chains
                .filter(line -> line.endsWith("1"))
                // filter for interesting binding modes with AMP/ATP or something like that present
                .filter(line -> Integer.parseInt(line.split(",")[11]) < 4)
                // skip bs structures
                .filter(line -> !line.contains("1jzs") // no AxP
                        && !line.contains("3vu8") // not mode 1-3
                        && !line.contains("4hwt") // no AxP
                        && !line.contains("5ah5") && !line.contains("3ziu") // naming fails
                        )
                .map(line -> line.split(","))
                .peek(split -> logger.info("fetching and parsing {}", split[2]))
                .collect(Collectors.toMap(split -> split[2] + "_" + split[13], split -> ProteinParser.source(split[2]).parse()));

        // extract key list
        List<String> ids = proteinMap.entrySet().stream()
                .map(Map.Entry::getKey)
                .collect(Collectors.toList());

        // extract ligands
        Map<String, AtomContainer> ligandMap = proteinMap.entrySet().stream()
                .collect(Collectors.toMap(Map.Entry::getKey, entry -> {
                    String id = entry.getKey();
                    logger.info("determining representative ligands {}", id);

                    Protein protein = entry.getValue();
                    GroupContainer representativeChain = protein.select()
                            .chainName(entry.getKey().split("_")[1])
                            .asGroupContainer();

                    AtomContainer ligand = protein.select()
                            .hetatms()
                            .asFilteredGroups()
                            // screen for ligands containing nitrogen
                            .filter(group -> group.getAtoms().stream()
                                    .map(Atom::getElement)
                                    .anyMatch(element -> element.equals(Element.N)))
                            .sorted(Comparator.comparingDouble((Group potentialLigand) -> computeMinimalDistanceToChain(potentialLigand, representativeChain)).reversed())
                            .findFirst()
                            .orElseThrow(() -> new NoSuchElementException("did not find any ligand for " + id));

                    logger.info("chosen ligand: {}", ligand);

                    return ligand;
                }));

        Set<String> atomNames = Stream
                .of("N1", "N3", "N7", "N9")
                .collect(Collectors.toSet());
        SVDSuperimposer svdSuperimposer = new SVDSuperimposer(atomNames, atomNames);

        // align with respect to ligand
        AtomContainer reference = ligandMap.get(ids.get(0));
        for(int alignmentIndex = 1; alignmentIndex < ids.size(); alignmentIndex++) {
            String id = ids.get(alignmentIndex);
            AtomContainer candidate = ligandMap.get(id);
            Protein protein = proteinMap.get(id);
            System.out.println(proteinMap.get(ids.get(0)) + " vs " + protein);
            StructureAlignmentResult alignmentResult = svdSuperimposer.align(reference, candidate);
            System.out.println(alignmentResult.getAlignmentScore());
            alignmentResult.transform(protein);
        }

        for (String id : ids) {
            String[] split = Files.lines(Paths.get(MAIN_TABLE_CATALYTIC_PATH))
                    .filter(line -> line.contains(id.split("_")[0].toLowerCase()))
                    .findFirst()
                    .get()
                    .split(",");

            GroupContainer chain = proteinMap.get(id).select()
                    .chainName(id.split("_")[1])
                    .aminoAcids()
                    .asGroupContainer();
            AtomContainer ligand = ligandMap.get(id);
            String output = chain.getPdbRepresentation() + System.lineSeparator() + ligand.getPdbRepresentation();
            Files.write(Paths.get("/home/bittrich/bs/" + split[0] + "-" + split[1] +  "-" + id + ".pdb"), output.getBytes());
        }
    }

    private static double computeMinimalDistanceToChain(Group potentialLigand, GroupContainer representativeChain) {
        double[] centroid = potentialLigand.calculate().centroid().getValue();

        return representativeChain.atoms()
                .map(Atom::calculate)
                .mapToDouble(coordinates -> coordinates.distanceFast(centroid))
                .min()
                .orElseThrow(NoSuchElementException::new);
    }
}
