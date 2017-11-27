package de.bioforscher.jstructure.contacts.collect.reconstruction.sathyapriya;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.contacts.ContactsConstants;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.stream.Collectors;

public class A05_ComputeTMAlignStatistics {
    public static void main(String[] args) {
        String output = ContactsConstants.lines(ContactsConstants.RECONSTRUCTION_DIRECTORY.resolve("ids.list"))
                .map(A05_ComputeTMAlignStatistics::handlePdbId)
                .collect(Collectors.joining(System.lineSeparator(),
                        "pdbId,sampling,mode,mapId,modelId,contacts,peraa,rmsd,tmscore" + System.lineSeparator(),
                        ""));

        ContactsConstants.write(ContactsConstants.RECONSTRUCTION_DIRECTORY.resolve("reconstruction-performance.csv"),
                output);
    }

    private static String handlePdbId(String pdbId) {
        Path target = ContactsConstants.RECONSTRUCTION_DIRECTORY.resolve("pdb").resolve(pdbId + ".pdb");
        return ContactsConstants.walk(ContactsConstants.RECONSTRUCTION_DIRECTORY.resolve("reconstructions"))
                .filter(path -> !Files.isDirectory(path))
                .filter(path -> path.toFile().getName().startsWith(pdbId) && path.toFile().getName().endsWith(".pdb") && path.toFile().getName().contains("_model"))
                .map(model -> handlePair(target, model))
                .collect(Collectors.joining(System.lineSeparator()));
    }

    private static String handlePair(Path target, Path model) {
        String modelName = model.toFile().getName().split("\\.")[0];
        System.out.println(modelName);

        String modelId = modelName.split("model")[1];
        String sampling = model.getParent().getParent().getParent().toFile().getName();
        String mode = model.getParent().getParent().toFile().getName().split("-")[1];
        String pdbId = model.toFile().getName().split("_")[0];

        String mapId;
        Path contactMapPath;
        if(sampling.equals("p100")) {
            mapId = "1";
            contactMapPath = ContactsConstants.RECONSTRUCTION_DIRECTORY.resolve("maps")
                    .resolve("p100")
                    .resolve(pdbId + "_A-" + mode + ".rr");
        } else {
            mapId = model.getParent().getParent().toFile().getName().split("-")[2];
            contactMapPath = ContactsConstants.RECONSTRUCTION_DIRECTORY.resolve("maps")
                    .resolve(sampling)
                    .resolve(pdbId + "_A-" + mode + "-" + mapId + ".rr");
        }

        int numberOfContacts = (int) ContactsConstants.lines(contactMapPath).count() - 1;
        int numberOfResidues = ContactsConstants.lines(ContactsConstants.RECONSTRUCTION_DIRECTORY.resolve("maps").resolve(pdbId + "_A.fasta"))
                .skip(1)
                .mapToInt(String::length)
                .findFirst()
                .getAsInt();

        double[] scores = computeAlignmentScores(target, model);
        return pdbId + "," +
                sampling + "," +
                mode + "," +
                mapId + "," +
                modelId + "," +
                numberOfContacts + "," +
                StandardFormat.format(numberOfContacts / (double) numberOfResidues) + "," +
                scores[0] + "," +
                scores[1];
    }

    private static double[] computeAlignmentScores(Path target, Path model) {
        try {
            ProcessBuilder processBuilder = new ProcessBuilder("tmalign",
                    target.toFile().getAbsolutePath(),
                    model.toFile().getAbsolutePath());

            Process process = processBuilder.start();
            List<String> lines = new BufferedReader(new InputStreamReader(process.getInputStream())).lines().collect(Collectors.toList());
            process.waitFor();

            double rmsd = 0;
            double tmscore = 0;

            for(String line : lines) {
                System.out.println(line);
                if(line.startsWith("Aligned length=")) {
                    rmsd = Double.valueOf(line.split("RMSD=")[1].trim().split(",")[0].trim());
                }
                if(line.startsWith("TM-score=")) {
                    tmscore = Double.valueOf(line.split("TM-score=")[1].trim().split("\\(")[0].trim());
                }
            }

            return new double[] { rmsd, tmscore };
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }
}
