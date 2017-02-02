package aars;

import de.bioforscher.jstructure.alignment.multiple.MultipleStructureAlignment;
import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceAreaCalculator;
import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProviderRegistry;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Element;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.family.AminoAcidFamily;
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

/**
 * Run MSA on aars data set.
 * Created by bittrich on 2/1/17.
 */
public class MultipleStructureAlignmentRunner {
    public static void main(String[] args) throws IOException {
        new MultipleStructureAlignmentRunner().performMultipleStructureAlignment();
    }

    private static final Logger logger = LoggerFactory.getLogger(MultipleStructureAlignmentRunner.class);

    private static final String homePath = System.getProperty("user.home");
    private static final String basePath = homePath + "/git/aars_analysis/data/";
    private static final AbstractFeatureProvider asaCalculator = FeatureProviderRegistry.resolve(AccessibleSurfaceAreaCalculator.ACCESSIBLE_SURFACE_AREA);
    private static final double ASA_CUTOFF = 0.16; // 10.1002/prot.340200303/epdf

    private void performMultipleStructureAlignment() throws IOException {
        int arginine1 = 884;
        int arginine2 = 1765;

        logger.info("loading structures...");
        List<GroupContainer> chains = Files.lines(Paths.get(basePath + "geometry/argtweezer_geometry.tsv"))
                .limit(10)
                .filter(line -> !line.startsWith("id"))
                .map(line -> line.split("\t"))
                .map(split -> split[0])
                .map(id -> {
                    logger.info("loading, parsing and computing asa for {}", id);

                    Protein protein = ProteinParser.parsePDBFile(basePath + "msa/C2/renumbered_structures/" + id.split("_")[0].toLowerCase() + "_renum.pdb");

                    Protein container = (Protein) Selection.on(protein)
                            .chainName(id.split("_")[1])
                            .asChainContainer();
                    container.setIdentifier(id);

                    asaCalculator.process(container);

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

                String aa = determineType(alignedFullStructure);

                Files.write(Paths.get(homePath + "/multiple-alignment/core/" + aa + "-" + numberOfAlignedGroups + "-" + alignedCore.getIdentifier() + ".pdb"), alignedCore.composePDBRecord().getBytes());
                Files.write(Paths.get(homePath + "/multiple-alignment/full/" + aa + "-" + alignedFullStructure.getIdentifier() + ".pdb"), alignedFullStructure.composePDBRecord().getBytes());
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
        groupsAround(chains, ampCentroid, homePath + "/multiple-alignment/amp/");
        System.out.println("around AA");
        groupsAround(chains, aaCentroid, homePath + "/multiple-alignment/aa/");
    }

    private String determineType(GroupContainer container) throws IOException {
        return Files.lines(Paths.get(basePath + "aars_main_table_catalytic.csv"))
                .filter(line -> line.contains(container.getIdentifier().toLowerCase().split("_")[0]))
                .findFirst()
                .orElseThrow(() -> new IllegalArgumentException("no value present for '" + container.getIdentifier() + "'"))
                .split(",")[1];
    }

    private void groupsAround(List<GroupContainer> chains, double[] centroid, String path) {
        final double bindingSiteExtension = 10.0;

        chains.forEach(chain -> {
            GroupContainer surroundingGroups = Selection.on(chain)
                    .aminoAcids()
                    .distance(centroid, bindingSiteExtension)
                    .asGroupContainer();

            surroundingGroups.getGroups().removeIf(group -> group.getFeatureAsDouble(AccessibleSurfaceAreaCalculator.ACCESSIBLE_SURFACE_AREA) / AminoAcidFamily.valueOfIgnoreCase(group.getThreeLetterCode()).orElse(AminoAcidFamily.UNKNOWN).getMaximalAccessibleSurfaceArea() < ASA_CUTOFF);

            try {
                String id = chain.getIdentifier();
                String aa = determineType(surroundingGroups);

                System.out.println(id + " : " + aa + " : " + surroundingGroups.groups()
                        .map(Group::getIdentifier)
                        .collect(Collectors.joining(", ")));

                Files.write(Paths.get(path + aa + "-" + id + ".pdb"), surroundingGroups.composePDBRecord().getBytes());
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
        });
    }
}
