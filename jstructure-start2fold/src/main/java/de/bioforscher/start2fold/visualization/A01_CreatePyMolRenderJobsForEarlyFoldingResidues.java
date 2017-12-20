package de.bioforscher.start2fold.visualization;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.start2fold.Start2FoldConstants;
import de.bioforscher.start2fold.model.Start2FoldResidueAnnotation;
import de.bioforscher.start2fold.parser.Start2FoldXmlParser;

import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

public class A01_CreatePyMolRenderJobsForEarlyFoldingResidues {
    public static void main(String[] args) throws IOException {
        String pymolCommand = Files.lines(Start2FoldConstants.PANCSA_LIST)
                .map(A01_CreatePyMolRenderJobsForEarlyFoldingResidues::composePyMolCommand)
                .collect(Collectors.joining(System.lineSeparator(),
                        // global settings
                        "bg_color white" + System.lineSeparator() +
                                "unset depth_cue" + System.lineSeparator() +
                                "set ray_opaque_background, off" + System.lineSeparator() +
                                "bg_color white" + System.lineSeparator() +
                                // create custom color palette
                                "set_color efr=[" + StandardFormat.format(23 / 255.0) + ", " + StandardFormat.format(111 / 255.0) + ", " + StandardFormat.format(193 / 255.0) + "]" + System.lineSeparator() +
                                "set ray_trace_mode, 3" + System.lineSeparator(),
                        ""));
        Files.write(Start2FoldConstants.PYMOL_DIRECTORY.resolve("efr.pml"),
                pymolCommand.getBytes());
    }

    private static String composePyMolCommand(String line) {
        String[] split = line.split(";");
        String entryId = split[0];
        String pdbId = split[1];

        List<Integer> experimentIds = Pattern.compile(",")
                .splitAsStream(split[2].replaceAll("\\[", "").replaceAll("]", ""))
                .map(Integer::valueOf)
                .collect(Collectors.toList());

        Structure structure = StructureParser.source(pdbId).parse();
        Chain chain = structure.chains().findFirst().get();

        Start2FoldXmlParser.parseSpecificExperiment(chain,
                Start2FoldConstants.XML_DIRECTORY.resolve(entryId + ".xml"),
                experimentIds);

        List<Integer> earlyFoldingResidues = chain.aminoAcids()
                .filter(aminoAcid -> aminoAcid.getFeature(Start2FoldResidueAnnotation.class).isEarly())
                .map(AminoAcid::getResidueIdentifier)
                .map(ResidueIdentifier::getResidueNumber)
                .collect(Collectors.toList());

        return "delete all" + System.lineSeparator() +
                "fetch " + pdbId + ", async=0" + System.lineSeparator() +
                // hide non-relevant stuff
                "hide everything" + System.lineSeparator() +
                "show cartoon, chain A" + System.lineSeparator() +
                // decolor everything
                "color grey80" + System.lineSeparator() +
                "zoom (chain A)" + System.lineSeparator() +
                earlyFoldingResidues.stream()
                        .map(res -> "color efr, resi " + res)
                        .collect(Collectors.joining(System.lineSeparator())) + System.lineSeparator() +
                "ray" + System.lineSeparator() +
                "png " + Start2FoldConstants.PYMOL_DIRECTORY.resolve(entryId + "-efr.png") + System.lineSeparator();
    }
}
