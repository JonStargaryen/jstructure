package de.bioforscher.jstructure.contacts.collect.reconstruction.sathyapriya;

import de.bioforscher.jstructure.contacts.ContactsConstants;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.stream.Collectors;

public class A05_CreateTMAlignScript {
    public static void main(String[] args) {
        String output = ContactsConstants.lines(ContactsConstants.RECONSTRUCTION_DIRECTORY.resolve("ids.list"))
                .map(A05_CreateTMAlignScript::handlePdbId)
                .collect(Collectors.joining(System.lineSeparator()));

        ContactsConstants.write(ContactsConstants.RECONSTRUCTION_DIRECTORY.resolve("tmalign.sh"),
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
        System.out.println(model);

        String sampling = model.getParent().getParent().getParent().toFile().getName();
        String mapName = model.getParent().getParent().toFile().getName();

        Path out = ContactsConstants.RECONSTRUCTION_DIRECTORY.resolve("tmalign")
                .resolve(sampling)
                .resolve(mapName)
                .resolve(model.toFile().getName().split("\\.")[0] + ".tmalign");

        ContactsConstants.ensureDirectoriesExist(out);

        return composeTMAlignCall(target, model, out);
    }

    private static String composeTMAlignCall(Path target, Path model, Path out) {
        return "echo " + model.toFile().getAbsolutePath() + System.lineSeparator() +
                "/home/bittrich/programs/tmalign/tmalign " + target.toFile().getAbsolutePath() + " " + model.toFile().getAbsolutePath() + " > " + out.toFile().getAbsolutePath();
    }
}
