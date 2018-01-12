package de.bioforscher.start2fold.visualization;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.start2fold.Start2FoldConstants;

import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Optional;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class A03_CreatePyMolRenderJobsForFunctionalResidues {
    public static void main(String[] args) throws IOException {
        String pymolCommand = Files.lines(Start2FoldConstants.PANCSA_LIST)
                .map(A03_CreatePyMolRenderJobsForFunctionalResidues::composePyMolCommand)
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.joining(System.lineSeparator(),
                        // global settings
                        "bg_color white" + System.lineSeparator() +
                                "unset depth_cue" + System.lineSeparator() +
                                "set ray_opaque_background, off" + System.lineSeparator() +
                                "bg_color white" + System.lineSeparator() +
                                // create custom color palette
                                "set_color efr=[" + StandardFormat.format(244 / 255.0) + ", " + StandardFormat.format(114 / 255.0) + ", " + StandardFormat.format(22 / 255.0) + "]" + System.lineSeparator() +
                                "set ray_trace_mode, 3" + System.lineSeparator(),
                        ""));
        Files.write(Start2FoldConstants.PYMOL_DIRECTORY.resolve("functional.pml"),
                pymolCommand.getBytes());
    }

    private static Optional<String> composePyMolCommand(String line) {
        String[] split = line.split(";");
        String entryId = split[0];
        String pdbId = split[1];

        List<Integer> functionalResidues = Pattern.compile(",")
                .splitAsStream(split[5].replaceAll("\\[", "").replaceAll("]", ""))
                .flatMap(A03_CreatePyMolRenderJobsForFunctionalResidues::valueOf)
                .collect(Collectors.toList());

        // skip entries missing an annotation of functional residues
        if(functionalResidues.isEmpty()) {
            return Optional.empty();
        }

        return Optional.of("delete all" + System.lineSeparator() +
                "fetch " + pdbId + ", async=0" + System.lineSeparator() +
                // hide non-relevant stuff
                "hide everything" + System.lineSeparator() +
                "show cartoon, chain A" + System.lineSeparator() +
                // decolor everything
                "color grey80" + System.lineSeparator() +
                "zoom (chain A)" + System.lineSeparator() +
                functionalResidues.stream()
                        .map(res -> "color efr, resi " + res)
                        .collect(Collectors.joining(System.lineSeparator())) + System.lineSeparator() +
                "ray" + System.lineSeparator() +
                "png " + Start2FoldConstants.PYMOL_DIRECTORY.resolve(entryId + "-functional.png") + System.lineSeparator());
    }

    private static Stream<Integer> valueOf(String string) {
        if(string.contains("-")) {
            String[] split = string.split("-");
            return IntStream.rangeClosed(Integer.valueOf(split[0]), Integer.valueOf(split[1]))
                    .boxed();
        } else {
            return Stream.of(Integer.valueOf(string));
        }
    }
}
