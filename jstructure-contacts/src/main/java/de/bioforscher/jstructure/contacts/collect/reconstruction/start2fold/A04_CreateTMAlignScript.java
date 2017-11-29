package de.bioforscher.jstructure.contacts.collect.reconstruction.start2fold;

import de.bioforscher.jstructure.contacts.ContactsConstants;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Optional;
import java.util.stream.Collectors;

public class A04_CreateTMAlignScript {
    public static void main(String[] args) {
        String output = ContactsConstants.lines(ContactsConstants.START2FOLD_DIRECTORY.resolve("ids.list"))
                .map(A04_CreateTMAlignScript::handlePdbId)
                .collect(Collectors.joining(System.lineSeparator()));

        ContactsConstants.write(ContactsConstants.START2FOLD_DIRECTORY.resolve("tmalign.sh"),
                output);
    }

    private static String handlePdbId(String pdbId) {
        Path target = ContactsConstants.START2FOLD_DIRECTORY.resolve("pdb").resolve(pdbId + ".pdb");
        return ContactsConstants.walk(ContactsConstants.START2FOLD_DIRECTORY.resolve("reconstructions"))
                .filter(path -> !Files.isDirectory(path))
                .filter(path -> path.toFile().getName().startsWith(pdbId) && path.toFile().getName().endsWith(".pdb") && path.toFile().getName().contains("_model"))
                .map(model -> handlePair(target, model))
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.joining(System.lineSeparator()));
    }

    private static Optional<String> handlePair(Path target, Path model) {
        String modelName = model.toFile().getName().split("\\.")[0];
        System.out.println(modelName);

        String reconstructionName = model.getParent().getParent().toFile().getName();

        Path out = ContactsConstants.START2FOLD_DIRECTORY.resolve("tmalign")
                .resolve(reconstructionName)
                .resolve(model.toFile().getName().split("\\.")[0] + ".tmalign");

        ContactsConstants.ensureDirectoriesExist(out);

        // tmalign output file already exists - do not do anything
        if(Files.exists(out)) {
            return Optional.empty();
        }

        return Optional.of(composeTMAlignCall(target, model, out));
    }

    private static String composeTMAlignCall(Path target, Path model, Path out) {
        return "echo " + model.toFile().getAbsolutePath() + System.lineSeparator() +
                "/home/bittrich/programs/tmalign/tmalign " + target.toFile().getAbsolutePath() + " " + model.toFile().getAbsolutePath() + " > " + out.toFile().getAbsolutePath();
    }
}
