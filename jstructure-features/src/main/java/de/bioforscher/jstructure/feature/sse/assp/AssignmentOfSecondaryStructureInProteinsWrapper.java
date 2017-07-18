package de.bioforscher.jstructure.feature.sse.assp;

import de.bioforscher.jstructure.model.feature.ComputationException;
import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.feature.sse.SecondaryStructureElement;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.selection.IntegerRange;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

/**
 * Use the CLI of ASSP.
 * Be careful: this implementation is not thread-safe per se as it uses the external executable and has a state in the
 * face that it writes temporary files to a working directory - when protein names equal, multiple threads could infer.
 * Created by bittrich on 7/5/17.
 */
@FeatureProvider(provides = { GenericSecondaryStructure.class }, priority = 50)
public class AssignmentOfSecondaryStructureInProteinsWrapper extends AbstractFeatureProvider {
    /**
     * Tied to a specific directory as files should be placed next to executable.
     */
    private static final Path ASSP_DIRECTORY = Paths.get("/usr/local/bin/assp/");

    @Override
    protected void processInternally(Structure protein) {
        try {
            protein.aminoAcids()
                    .forEach(this::assignNeutralState);

            // create temp file next to executable as expected by ASSP
            String pdbId = protein.getProteinIdentifier().getFullName();
            String filename = pdbId + ".pdb";
            Path tmpPath = ASSP_DIRECTORY.resolve(filename);
            Files.write(tmpPath, protein.getPdbRepresentation().getBytes());

            ProcessBuilder processBuilder = new ProcessBuilder("./ASSP_lin.exe", "-i", filename);
            processBuilder.directory(ASSP_DIRECTORY.toFile());
            processBuilder.start().waitFor();

            Path outputPath = ASSP_DIRECTORY.resolve(pdbId + "_assp.out");
            while(!Files.exists(outputPath)) {
                System.out.println("waiting");
                Thread.sleep(10);
            }
            Files.lines(outputPath)
                    .forEach(line -> handleLine(protein, line));

            // clear all tmp files
            Files.delete(tmpPath);
            Files.delete(outputPath);
            Files.delete(ASSP_DIRECTORY.resolve(pdbId+ "_cont.out"));
            Files.delete(ASSP_DIRECTORY.resolve(pdbId + "_diff.out"));
            Files.delete(ASSP_DIRECTORY.resolve(pdbId + "_nh.out"));
        } catch (Exception e) {
            throw new ComputationException(e);
        }
    }

    private void assignNeutralState(AminoAcid aminoAcid) {
        aminoAcid.getFeatureContainer().addFeature(new GenericSecondaryStructure(this, SecondaryStructureElement.COIL));
    }

    private void handleLine(Structure protein, String line) {
        line = line.replaceAll("\\s+", " ");
        String[] split = line.split(" ");

        SecondaryStructureElement secondaryStructure = mapToSecondaryStructureElement(split[0]);

        Chain chain = protein.select()
                .chainName(split[3])
                .asChain();
        chain.select()
                .residueNumber(new IntegerRange(Integer.valueOf(split[4]), Integer.valueOf(split[7])))
                .asFilteredGroups()
                .map(Group::getFeatureContainer)
                .map(featureContainer -> featureContainer.getFeature(GenericSecondaryStructure.class))
                .forEach(ss -> ss.setSecondaryStructure(secondaryStructure));
    }

    private SecondaryStructureElement mapToSecondaryStructureElement(String string) {
        switch(string) {
            case "AlphaHelix":
                return SecondaryStructureElement.ALPHA_HELIX;
            case "ThreeHelix":
                return SecondaryStructureElement.THREE_TEN_HELIX;
            case "PiHelix":
                return SecondaryStructureElement.PI_HELIX;
            default:
                throw new IllegalArgumentException("unknown structure: " + string);
        }
    }
}
