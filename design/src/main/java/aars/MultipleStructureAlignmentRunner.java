package aars;

import de.bioforscher.jstructure.alignment.multiple.MultipleStructureAlignment;
import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceAreaCalculator;
import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProviderRegistry;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Element;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import static aars.AARSConstants.*;

/**
 * Run MSA on aars data set.
 * Created by bittrich on 2/1/17.
 */
public class MultipleStructureAlignmentRunner {
    public static void main(String[] args) throws IOException {
        new MultipleStructureAlignmentRunner().performMultipleStructureAlignment();
    }

    private static final Logger logger = LoggerFactory.getLogger(MultipleStructureAlignmentRunner.class);
    private static final AbstractFeatureProvider asaCalculator = FeatureProviderRegistry.resolve(AccessibleSurfaceAreaCalculator.ACCESSIBLE_SURFACE_AREA);
//    private static final double ASA_CUTOFF = 0.16; // 10.1002/prot.340200303/epdf
    private static final String RELATIVE_ACCESSIBLE_SURFACE_AREA = "RELATIVE_ACCESSIBLE_SURFACE_AREA";

    private void performMultipleStructureAlignment() throws IOException {
        int arginine1 = 884;
        int arginine2 = 1765;

        logger.info("loading structures...");
        List<GroupContainer> chains = Files.lines(Paths.get(ARGTWEEZER_GEOMETRY_PATH))
                .filter(line -> !line.startsWith("id"))
                .map(line -> line.split("\t"))
                .map(split -> split[0])
                .map(id -> {
                    logger.info("loading, parsing and computing asa for {}", id);

                    Protein protein = ProteinParser.parsePDBFile(RENUMBERED_STRUCTURES_C2_PATH + id.split("_")[0].toLowerCase() + "_renum.pdb");

                    Protein container = (Protein) Selection.on(protein)
                            .chainName(id.split("_")[1])
                            .asChainContainer();
                    container.setIdentifier(id);

//                    asaCalculator.process(container);
//
//                    // assign each rasa value as bfactor
//                    container.aminoAcids().forEach(group -> {
//                        double rasa = group.getFeatureAsDouble(AccessibleSurfaceAreaCalculator.ACCESSIBLE_SURFACE_AREA) / AminoAcidFamily.valueOfIgnoreCase(group.getThreeLetterCode()).orElse(AminoAcidFamily.UNKNOWN).getMaximalAccessibleSurfaceArea();
//                        group.setFeature(RELATIVE_ACCESSIBLE_SURFACE_AREA, rasa);
//                        group.atoms().forEach(atom -> atom.setBfactor((float) rasa));
//                    });

                    return container;
                })
                .collect(Collectors.toList());

        MultipleStructureAlignment multipleStructureAlignment = new MultipleStructureAlignment();
        multipleStructureAlignment.align(chains, arginine1, arginine2);
        Map<GroupContainer, GroupContainer> alignedContainers = multipleStructureAlignment.getAlignedContainerMap();

        // compute distance matrix for MDS visualization
//        List<String> identifiers = alignedContainers.entrySet().stream()
//                .map(Map.Entry::getKey)
//                .map(GroupContainer::getIdentifier)
//                .filter(identifier -> {
//                    try {
//                        return Files.lines(Paths.get("/home/bittrich/git/aars_analysis/data/aars_main_table_catalytic.csv"))
//                                .filter(line -> line.contains(identifier.split("_")[0].toLowerCase()))
//                                .anyMatch(line -> line.endsWith("1"));
//                    } catch (IOException e) {
//                        throw new UncheckedIOException(e);
//                    }
//                })
//                .collect(Collectors.toList());
//        List<String> outputLines = new ArrayList<>();
//        for(int rowIndex = 0; rowIndex < identifiers.size(); rowIndex++) {
//            List<String> rmsds = new ArrayList<>();
//            for (String identifier : identifiers) {
//                double rmsd;
//                Pair<String, String> key = new Pair<>(identifiers.get(rowIndex), identifier);
//                if (identifiers.get(rowIndex).equals(identifier)) {
//                    rmsd = 0;
//                } else {
//                    rmsd = LinearAlgebraAtom.calculateRmsd(alignedContainers.entrySet().stream()
//                                    .filter(pair -> pair.getKey().getIdentifier().equals(key.getLeft()))
//                                    .findFirst()
//                                    .get()
//                                    .getValue(),
//                            alignedContainers.entrySet().stream()
//                                    .filter(pair -> pair.getKey().getIdentifier().equals(key.getRight()))
//                                    .findFirst()
//                                    .get()
//                                    .getValue());
//                }
//                rmsds.add(String.valueOf(rmsd));
//            }
//            String outputLine = rmsds.stream().collect(Collectors.joining(",", identifiers.get(rowIndex) + ",", ""));
//            outputLines.add(outputLine);
//        }
//        String output = outputLines.stream()
//                .collect(Collectors.joining(System.lineSeparator(), identifiers.stream()
//                        .collect(Collectors.joining(",", ",", System.lineSeparator())), System.lineSeparator()));
//        Files.write(Paths.get("/home/bittrich/bindingsite-rmsd.csv"), output.getBytes());
//        // update labels.txt
//        String newLabels = Files.lines(Paths.get("/home/bittrich/bindingsite-rmsd.csv"))
//                .filter(line -> !line.startsWith(","))
//                .map(line -> line.split(",")[0])
//                .map(id -> id.split("_")[0])
//                .map(String::toLowerCase)
//                .map(id -> {
//                    try {
//                        return Files.lines(Paths.get("/home/bittrich/git/aars_analysis/data/aars_main_table_catalytic.csv"))
//                                .filter(line -> line.contains(id))
//                                .findFirst()
//                                .get()
//                                .split(",")[1];
//                    } catch (IOException e) {
//                        throw new UncheckedIOException(e);
//                    }
//                })
//                .map(aa -> {
//                    System.out.println(aa);
//                    try {
//                        return Files.lines(Paths.get("/home/bittrich/git/aars_analysis/data/pocket_align/labels.txt"))
//                                .filter(line -> line.startsWith(aa))
//                                .findFirst()
//                                .get();
//                    } catch (IOException e) {
//                        throw new UncheckedIOException(e);
//                    }
//                })
//                .collect(Collectors.joining(System.lineSeparator(), "label,color", System.lineSeparator()));
//        Files.write(Paths.get("/home/bittrich/git/aars_analysis/data/pocket_align/labels_binding_site.txt"), newLabels.getBytes());
//        System.exit(0);
        // end: compute distance matrix for MDS visualization

        alignedContainers.entrySet().forEach(entry -> {
            try {
                GroupContainer alignedFullStructure = entry.getKey();
                GroupContainer alignedCore = entry.getValue();
                int numberOfAlignedGroups = alignedCore.getGroups().size();

                String aa = determineType(alignedFullStructure);

                Files.createDirectories(Paths.get(HOME_PATH + "/multiple-alignment/core/"));
                Files.write(Paths.get(HOME_PATH + "/multiple-alignment/core/" + aa + "_" + numberOfAlignedGroups + "-" + alignedCore.getIdentifier() + ".pdb"), alignedCore.composePDBRecord().getBytes());
                Files.createDirectories(Paths.get(HOME_PATH + "/multiple-alignment/full/"));
                Files.write(Paths.get(HOME_PATH + "/multiple-alignment/full/" + aa + "_" + alignedFullStructure.getIdentifier() + ".pdb"), alignedFullStructure.composePDBRecord().getBytes());
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
        });

        // compute centroids in order to extract both binding sites
        List<double[]> ampCoordinates = alignedContainers.entrySet().stream()
                .map(Map.Entry::getKey)
                .flatMap(container -> Selection.on(container)
                        .hetatms()
                        .groupName("AMP")
                        // remove some nitrogen far away
                        .negationModeEnter()
                        .water()
                        .pdbSerial(6373, 6375, 6378, 6379, 6381)
                        .negationModeLeave()
                        .element(Element.N)
                        .asFilteredAtoms())
                .map(Atom::getCoordinates)
                .collect(Collectors.toList());
        double[] ampCentroid = LinearAlgebra3D.divide(ampCoordinates.stream()
                .reduce(new double[3], LinearAlgebra3D::add, LinearAlgebra3D::add), ampCoordinates.size());
        logger.info("AMP-centroid: " + Arrays.toString(ampCentroid));

        List<double[]> aaCoordinates = alignedContainers.entrySet().stream()
                .map(Map.Entry::getKey)
                .flatMap(container -> Selection.on(container)
                        .groupName("A5A", "AMO", "GAP", "HAM", "HSS", "KAA", "FA5", "5CA", "P5A", "YLY", "N0B", "YLA", "TSB")
                        .atomName("CA")
                        .asFilteredAtoms())
                .map(Atom::getCoordinates)
                .collect(Collectors.toList());
        double[] aaCentroid = LinearAlgebra3D.divide(aaCoordinates.stream()
                .reduce(new double[3], LinearAlgebra3D::add, LinearAlgebra3D::add), aaCoordinates.size());
        logger.info("AA-centroid: " + Arrays.toString(aaCentroid));

        System.out.println("around AMP");
        groupsAround(chains, ampCentroid, HOME_PATH + "/multiple-alignment/amp/", 0);
        groupsAround(chains, ampCentroid, HOME_PATH + "/multiple-alignment/amp-rasa-16/",  0.16);
//        groupsAround(chains, ampCentroid, HOME_PATH + "/multiple-alignment/amp-rasa-36/", 0.36);
//        groupsAround(chains, ampCentroid, HOME_PATH + "/multiple-alignment/amp-rasa-56/", 0.56);
//        groupsAround(chains, ampCentroid, HOME_PATH + "/multiple-alignment/amp-rasa-76/", 0.76);
        System.out.println("around AA");
        groupsAround(chains, aaCentroid, HOME_PATH + "/multiple-alignment/aa/",  0);
        groupsAround(chains, aaCentroid, HOME_PATH + "/multiple-alignment/aa-rasa-16/", 0.16);
//        groupsAround(chains, ampCentroid, HOME_PATH + "/multiple-alignment/aa-rasa-36/", 0.36);
//        groupsAround(chains, ampCentroid, HOME_PATH + "/multiple-alignment/aa-rasa-56/", 0.56);
//        groupsAround(chains, ampCentroid, HOME_PATH + "/multiple-alignment/aa-rasa-76/", 0.76);
    }

    private String determineType(GroupContainer container) throws IOException {
        return Files.lines(Paths.get(MAIN_TABLE_CATALYTIC_PATH))
                .filter(line -> line.contains(container.getIdentifier().toLowerCase().split("_")[0]))
                .findFirst()
                .orElseThrow(() -> new IllegalArgumentException("no value present for '" + container.getIdentifier() + "'"))
                .split(",")[1];
    }

    private void groupsAround(List<GroupContainer> chains, double[] centroid, String path, double minimalRasa) {
        final double bindingSiteExtension = 10.0;

        chains.forEach(chain -> {
            GroupContainer surroundingGroups = Selection.on(chain)
                    .aminoAcids()
                    .distance(centroid, bindingSiteExtension)
                    .cloneElements()
                    .asGroupContainer();

            if(minimalRasa > 0.0) {
                surroundingGroups.getGroups().removeIf(group -> group.getFeatureAsDouble(RELATIVE_ACCESSIBLE_SURFACE_AREA) < minimalRasa);
            }

            try {
                String id = chain.getIdentifier();
                String aa = determineType(surroundingGroups);
                Files.createDirectories(Paths.get(path));
                Files.write(Paths.get(path + aa + "_" + id + ".pdb"), surroundingGroups.composePDBRecord().getBytes());
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
        });
    }
}
