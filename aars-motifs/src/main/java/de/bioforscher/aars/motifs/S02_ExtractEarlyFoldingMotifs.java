package de.bioforscher.aars.motifs;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureCollectors;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.testutil.FileUtils;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

public class S02_ExtractEarlyFoldingMotifs {
    private static final double CUTOFF = 0.163;

    public static void main(String[] args) {
        FileUtils.list(Paths.get("/home/bittrich/git/phd_sb_repo/data/aars/EFoldMine_code/pdb-motifs/ef-predictions/"))
                .parallel()
                .forEach(earlyFoldingPredictionPath -> {
                    String id = earlyFoldingPredictionPath.toFile().getName().split("\\.")[0];
                    System.out.println(id);
                    Path motifPath = earlyFoldingPredictionPath.getParent().getParent().resolve("motifs").resolve(id + ".pdb");

                    String[] split = id.split("_");
                    String pdbId = split[0];
                    String chainId = split[1];

                    try {
                        Structure structure = StructureParser.fromPdbId(pdbId).parse();
                        Chain chain = structure.select().chainId(chainId).asChain();
                        List<AminoAcid> aminoAcids = chain.getAminoAcids();
                        List<String> earlyFoldingPredictionLines = Files.readAllLines(earlyFoldingPredictionPath);
                        List<Double> earlyFoldingPredictions = earlyFoldingPredictionLines.stream()
                                .skip(1)
                                .filter(line -> !line.startsWith("*"))
                                .map(line -> line.split("\\s+")[1])
                                .map(Double::valueOf)
                                .collect(Collectors.toList());

                        if (aminoAcids.size() != earlyFoldingPredictions.size()) {
                            System.err.println("erroneous data for " + id);
                            return;
                        }

                        StringJoiner output = new StringJoiner(System.lineSeparator());
                        output.add(structure.getHeader());

                        List<AminoAcid> earlyFoldingResidues = new ArrayList<>();
                        for (int i = 0; i < aminoAcids.size(); i++) {
                            if (earlyFoldingPredictions.get(i) > CUTOFF) {
                                continue;
                            }

                            // is predicted as early folding
                            earlyFoldingResidues.add(aminoAcids.get(i));
                        }

                        Structure earlyFoldingStructure = earlyFoldingResidues.stream().collect(StructureCollectors.toIsolatedStructure());
                        FileUtils.write(motifPath, earlyFoldingStructure.getPdbRepresentation());
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                });
    }
}
