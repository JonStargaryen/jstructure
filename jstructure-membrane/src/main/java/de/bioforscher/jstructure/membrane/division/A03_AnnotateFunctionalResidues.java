package de.bioforscher.jstructure.membrane.division;

import de.bioforscher.jstructure.feature.mapping.SiftsMappingAnnotator;
import de.bioforscher.jstructure.feature.topology.MembraneContainer;
import de.bioforscher.jstructure.feature.topology.OrientationsOfProteinsInMembranesAnnotator;
import de.bioforscher.jstructure.membrane.MembraneConstants;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;

import java.nio.file.Path;
import java.util.Optional;
import java.util.stream.Collectors;

public class A03_AnnotateFunctionalResidues {
    private static final UniProtAnnotator UNI_PROT_ANNOTATOR = new UniProtAnnotator();
    private static final OrientationsOfProteinsInMembranesAnnotator ORIENTATIONS_OF_PROTEINS_IN_MEMBRANES_ANNOTATOR = new OrientationsOfProteinsInMembranesAnnotator();

    public static void main(String[] args) {
        String output = MembraneConstants.list(MembraneConstants.DIVISION_DIRECTORY.resolve("dynamine"))
                .limit(4)
                .map(A03_AnnotateFunctionalResidues::handleFile)
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.joining(System.lineSeparator()));

        MembraneConstants.write(MembraneConstants.DIVISION_DIRECTORY.resolve("functional.csv"),
                output);
    }

    private static Optional<String> handleFile(Path path) {
        try {
            String[] split = path.toFile().getName().split("_");
            String pdbId = split[0];
            String chainId = split[1];
            String id = pdbId + "_" + chainId;
            System.out.println(id);

            String uniProtId = SiftsMappingAnnotator.getLinesForPdbId(pdbId)
                    .filter(line -> line[1].equals(chainId))
                    .findFirst()
                    .get()[2];
            System.out.println(uniProtId);

            Structure structure = StructureParser.source(pdbId)
                    .minimalParsing(true)
                    .parse();
            Chain chain = structure.select()
                    .chainName(chainId)
                    .asChain();

            UniProtAnnotationContainer container = UNI_PROT_ANNOTATOR.process(uniProtId);

            System.out.println(container.getActiveSites());
//            System.out.println(container.getTransmembraneRegions());

            ORIENTATIONS_OF_PROTEINS_IN_MEMBRANES_ANNOTATOR.process(structure);
            chain.getFeature(MembraneContainer.class)
                    .getTransMembraneSubunits()
                    .stream()
                    .forEach(System.out::println);

            return Optional.of(id + ";" +
                    uniProtId);
        } catch (Exception e) {
            return Optional.empty();
        }
    }
}
