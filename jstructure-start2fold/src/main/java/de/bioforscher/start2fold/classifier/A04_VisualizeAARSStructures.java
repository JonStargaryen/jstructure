package de.bioforscher.start2fold.classifier;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.mathematics.LinearAlgebra;
import de.bioforscher.jstructure.mathematics.SetOperations;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.start2fold.Start2FoldConstants;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class A04_VisualizeAARSStructures {
    public static void main(String[] args) throws IOException {
        String pymolCommand = Files.lines(Paths.get("/home/bittrich/git/aars_data/dataset.csv"))
                .filter(line -> !line.startsWith("identifier"))
                .map(A04_VisualizeAARSStructures::composePyMolCommand)
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.joining(System.lineSeparator(),
                        // global settings
                        "bg_color white" + System.lineSeparator() +
                                "unset depth_cue" + System.lineSeparator() +
                                "set ray_opaque_background, off" + System.lineSeparator() +
                                "bg_color white" + System.lineSeparator() +
                                // create custom color palette
                                "set_color func=[" + StandardFormat.format(244 / 255.0) + ", " + StandardFormat.format(114 / 255.0) + ", " + StandardFormat.format(22 / 255.0) + "]" + System.lineSeparator() +
                                "set_color efr=[" + StandardFormat.format(23 / 255.0) + ", " + StandardFormat.format(111 / 255.0) + ", " + StandardFormat.format(193 / 255.0) + "]" + System.lineSeparator() +
                                "set ray_trace_mode, 3" + System.lineSeparator(),
                        ""));
        Files.write(Start2FoldConstants.DATA_DIRECTORY.resolve("classifier").resolve("aars").resolve("aars.pml"),
                pymolCommand.getBytes());
    }

    private static Optional<String> composePyMolCommand(String line) {
        String[] split = line.split(",");
        String id = split[0];
        String pdbId = id.split("_")[0];
        String chainId = id.split("_")[1];
        String clazz = split[1];

        // skip non-representative structures
        List<String> fastaLines;
        try(Stream<String> lines = Files.lines(Paths.get("/home/bittrich/git/aars_data/T04_representative_sequences/C" + clazz + "_representatives_cluster.fasta"))) {
            fastaLines = lines.collect(Collectors.toList());
        } catch (IOException e) {
//            e.printStackTrace();
            return Optional.empty();
        }

        if(fastaLines.stream().noneMatch(l -> l.equals(">" + pdbId + "_" + chainId))) {
            return Optional.empty();
        }

//        System.out.println(line);

        Structure originalStructure = StructureParser.fromPath(Paths.get("/home/bittrich/git/aars_data/T06_renumbered_structures/C" + clazz + "/" + pdbId + "_" + chainId + ".pdb"))
                .parse();
        Chain originalChain = originalStructure.select()
                .chainName(chainId)
                .asChain();
        Structure renumberedStructure = StructureParser.fromPath(Paths.get("/home/bittrich/git/aars_data/T06_renumbered_structures/C" + clazz + "/renum/" + pdbId + "_" + chainId + "_renum.pdb"))
                .parse();
        Chain renumberedChain = renumberedStructure.select()
                .chainName(chainId)
                .asChain();

        List<String> earlyLines;
        try(Stream<String> lines = Files.lines(Start2FoldConstants.DATA_DIRECTORY.resolve("classifier")
                .resolve("aars")
                .resolve("out")
                .resolve(pdbId + ".out"))) {
            earlyLines = lines.collect(Collectors.toList());
        } catch (IOException e) {
//            e.printStackTrace();
            return Optional.empty();
        }
        List<Integer> earlyFoldingResidues = earlyLines.stream()
                // ignore some header lines
                .filter(out -> !out.startsWith("chain"))
                // explicitly search for EFR
                .filter(out -> out.endsWith("early"))
                // search for correct chain
                .filter(out -> out.startsWith(chainId))
                .map(out -> out.split(","))
                // ignore most potential header lines
                .filter(outSplit -> outSplit.length > 0)
                // extract residue numbers
                .map(outSplit -> outSplit[1])
                .mapToInt(Integer::valueOf)
                .boxed()
                .collect(Collectors.toList());

        List<String> functionalLines;
        try(Stream<String> lines = Files.lines(Paths.get("/home/bittrich/git/aars_data/T09_interactions/C" + clazz + "/contacts_per_structure.txt"))) {
            functionalLines = lines.collect(Collectors.toList());
        } catch (IOException e) {
//            e.printStackTrace();
            return Optional.empty();
        }
        List<Integer> functionalResidues = functionalLines.stream()
                .filter(l -> l.startsWith(pdbId + "_" + chainId))
                .map(l -> l.split("\\[")[1].split("]")[0])
                .map(l -> l.split(","))
                .flatMap(Stream::of)
                .mapToInt(Integer::valueOf)
                .mapToObj(renumberedResidueNumber -> mapToOriginalResidueNumber(renumberedResidueNumber, originalChain, renumberedChain))
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.toList());

        if(earlyFoldingResidues.isEmpty() || functionalResidues.isEmpty()) {
            return Optional.empty();
        }

        System.out.println("EFR: " + earlyFoldingResidues);
        System.out.println("func: " + functionalResidues);
        System.out.println("overlap: " + SetOperations.createIntersectionSet(earlyFoldingResidues, functionalResidues));


        return Optional.of("delete all" + System.lineSeparator() +
                "fetch " + pdbId + ", async=0" + System.lineSeparator() +
                // hide non-relevant stuff
                "hide everything" + System.lineSeparator() +
                "show cartoon, chain " + chainId + System.lineSeparator() +
                // decolor everything
                "color grey80" + System.lineSeparator() +
                "zoom (chain " + chainId + ")" + System.lineSeparator() +
                earlyFoldingResidues.stream()
                        .map(res -> "color efr, resi " + res)
                        .collect(Collectors.joining(System.lineSeparator())) + System.lineSeparator() +
                functionalResidues.stream()
                        .map(res -> "color func, resi " + res)
                        .collect(Collectors.joining(System.lineSeparator())) + System.lineSeparator() +
                "ray" + System.lineSeparator() +
                "png " + Start2FoldConstants.DATA_DIRECTORY.resolve("classifier").resolve("aars").resolve("png").resolve(id + ".png") + System.lineSeparator());
    }

    private static Optional<Integer> mapToOriginalResidueNumber(int renumberedResidueNumber, Chain originalChain, Chain renumberedChain) {
        try {
            AminoAcid renumberedAminoAcid = renumberedChain.select()
                    .residueNumber(renumberedResidueNumber)
                    .asAminoAcid();
            LinearAlgebra.PrimitiveDoubleArrayLinearAlgebra renumberedPosition = renumberedAminoAcid.calculate()
                    .centroid();

            return originalChain.aminoAcids()
                    .filter(originalAminoAcid -> originalAminoAcid.calculate()
                            .centroid()
                            .distance(renumberedPosition) < 0.01)
                    .findFirst()
                    .map(Group::getResidueIdentifier)
                    .map(ResidueIdentifier::getResidueNumber);
        } catch (Exception e) {
            return Optional.empty();
        }
    }
}
