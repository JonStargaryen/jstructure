package aars;

import de.bioforscher.jstructure.alignment.AlignmentResult;
import de.bioforscher.jstructure.alignment.SVDSuperimposer;
import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Element;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import static aars.AARSConstants.MAIN_TABLE_CATALYTIC_PATH;

/**
 *
 * Aligns binding sites with respect to a present ligand.
 * Created by bittrich on 2/7/17.
 */
public class BindingSiteAligner {
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
                .filter(line -> !line.contains("1jzs"))
                .map(line -> line.split(","))
                .peek(split -> logger.info("fetching and parsing {}", split[2]))
                .collect(Collectors.toMap(split -> split[2] + "_" + split[13], split -> ProteinParser.parseProteinById(split[2])));

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
                    GroupContainer representativeChain = Selection.on(protein)
                            .chainName(entry.getKey().split("_")[1])
                            .asGroupContainer();

                    AtomContainer ligand = Selection.on(protein)
                            .hetatms()
                            .asFilteredGroups()
                            // screen for ligands containing nitrogen
                            .filter(group -> group.getAtoms().stream()
                                    .map(Atom::getElement)
                                    .anyMatch(element -> element.equals(Element.N)))
                            .sorted(Comparator.comparingDouble((Group potentialLigand) -> computeMimimalDistanceToChain(potentialLigand, representativeChain)).reversed())
                            .findFirst()
                            .orElseThrow(() -> new NoSuchElementException("did not find any ligand for " + id));

                    logger.info("chosen ligand: {}", ligand);

                    return ligand;
                }));

        SVDSuperimposer svdSuperimposer = new SVDSuperimposer(Collections.emptySet(), Stream.of("N", "C", "O")
                .flatMap(element -> IntStream.range(0, 11).mapToObj(i -> element + i)).collect(Collectors.toSet()));

        // align with respect to ligand
        AtomContainer reference = ligandMap.get(ids.get(0));
        for(int alignmentIndex = 1; alignmentIndex < ids.size(); alignmentIndex++) {
            String id = ids.get(alignmentIndex);
            AtomContainer candidate = ligandMap.get(id);
            AlignmentResult alignmentResult = svdSuperimposer.align(reference, candidate);
            System.out.println(reference + " vs " + candidate + " : " + alignmentResult.getAlignmentScore());
            Protein protein = proteinMap.get(id);
            alignmentResult.transform(protein);
        }

        for(int alignmentIndex = 0; alignmentIndex < ids.size(); alignmentIndex++) {
            String id = ids.get(alignmentIndex);
            Files.write(Paths.get("/home/bittrich/bs/" + id + ".pdb"), proteinMap.get(id).composePDBRecord().getBytes());
        }
    }

    private static double computeMimimalDistanceToChain(Group potentialLigand, GroupContainer representativeChain) {
        double[] centroid = potentialLigand.atoms()
                .map(Atom::getCoordinates)
                .reduce(new double[3], LinearAlgebra3D::add, LinearAlgebra3D::add);

        return representativeChain.atoms()
                .map(Atom::getCoordinates)
                .mapToDouble(coordinates -> LinearAlgebra3D.distanceFast(centroid, coordinates))
                .min()
                .orElseThrow(NoSuchElementException::new);
    }
}
