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
import org.jsoup.Jsoup;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

public class A02_CreatePyMolRenderJobsForStrongResidues {
    public static void main(String[] args) throws IOException {
        String pymolCommand = Files.list(Start2FoldConstants.XML_DIRECTORY)
                // ignore the file whose sequences cannot be aligned
                .filter(path -> !path.toFile().getName().contains("STF0034"))
                .map(A02_CreatePyMolRenderJobsForStrongResidues::composePyMolCommand)
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
        Files.write(Start2FoldConstants.PYMOL_DIRECTORY.resolve("strong.pml"),
                pymolCommand.getBytes());
    }

    private static Optional<String> composePyMolCommand(Path path) {
        try {
            String entryId = path.toFile().getName().split("\\.")[0];
            String pdbId = Jsoup.parse(path.toFile(), "UTF-8").getElementsByTag("protein").attr("pdb_id");

            Structure structure = StructureParser.source(pdbId).parse();
            Chain chain = structure.chains().findFirst().get();

            Start2FoldXmlParser.parse(chain,
                    Start2FoldConstants.XML_DIRECTORY.resolve(entryId + ".xml"));

            List<Integer> strongResidues = chain.aminoAcids()
                    .filter(aminoAcid -> aminoAcid.getFeature(Start2FoldResidueAnnotation.class).isStrong())
                    .map(AminoAcid::getResidueIdentifier)
                    .map(ResidueIdentifier::getResidueNumber)
                    .collect(Collectors.toList());

            if (strongResidues.isEmpty()) {
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
                    strongResidues.stream()
                            .map(res -> "color efr, resi " + res)
                            .collect(Collectors.joining(System.lineSeparator())) + System.lineSeparator() +
                    "ray" + System.lineSeparator() +
                    "png " + Start2FoldConstants.PYMOL_DIRECTORY.resolve(entryId + "-strong.png") + System.lineSeparator());
        } catch (IOException e) {
            return Optional.empty();
        }
    }
}
