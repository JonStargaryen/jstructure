package de.bioforscher.start2fold.collect;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.start2fold.Start2FoldConstants;
import de.bioforscher.start2fold.model.Start2FoldResidueAnnotation;
import de.bioforscher.start2fold.parser.EvolutionaryCouplingParser;
import de.bioforscher.start2fold.parser.Start2FoldXmlParser;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.select.Elements;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * Assess whether couplings tend to be local or long-range contacts.
 *
 * Conclusion: all couplings are by definition long-range contacts
 */
public class A07_WriteCouplingRangeCsv {
    public static void main(String[] args) throws IOException {
        String output = Files.lines(Start2FoldConstants.PANCSA_LIST)
                .map(A07_WriteCouplingRangeCsv::handleLine)
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.joining(System.lineSeparator(),
                        "pdb,chain,res,aa,folds,range,plm_score" + System.lineSeparator(),
                        ""));

        Start2FoldConstants.write(Start2FoldConstants.STATISTICS_DIRECTORY.resolve("couplings.csv"),
                output);
    }

    private static Optional<String> handleLine(String line) {
        try {
            System.out.println(line);
            String[] split = line.split(";");
            String entryId = split[0];
            String pdbId = split[1];
            List<Integer> experimentIds = Pattern.compile(",")
                    .splitAsStream(split[2].replaceAll("\\[", "").replaceAll("]", ""))
                    .map(Integer::valueOf)
                    .collect(Collectors.toList());

            Structure structure = StructureParser.fromPdbId(pdbId).parse();
            Chain chain = structure.chains().findFirst().get();

            Start2FoldXmlParser.parseSpecificExperiment(chain,
                    Start2FoldConstants.XML_DIRECTORY.resolve(entryId + ".xml"),
                    experimentIds);

            EvolutionaryCouplingParser.parseHotSpotFile(chain,
                    Start2FoldConstants.COUPLING_DIRECTORY.resolve(entryId.toUpperCase() + "_hs.html"));

            List<AminoAcid> earlyFoldingResidues = chain.aminoAcids()
                    .filter(aminoAcid -> aminoAcid.getFeature(Start2FoldResidueAnnotation.class).isEarly())
                    .collect(Collectors.toList());

            Map<Integer, List<Double>> localPlmScores = new HashMap<>();
            Map<Integer, List<Double>> longRangePlmScores = new HashMap<>();

            Document hotSpotDocument = Jsoup.parse(Files.readAllLines(Paths.get("/home/bittrich/git/phd_sb_repo/data/start2fold/coupling/" + entryId + "_ec.html"))
                    .stream()
                    .collect(Collectors.joining(System.lineSeparator())));
            List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());

            for(int i = 0; i < aminoAcids.size(); i++) {
                localPlmScores.put(i, new ArrayList<>());
                longRangePlmScores.put(i, new ArrayList<>());
            }

            hotSpotDocument.getElementsByTag("tr").stream()
                    .skip(1)
                    .forEach(element -> {
                        Elements tds = element.getElementsByTag("td");
                        int residueNumber1 = Integer.valueOf(tds.get(2).text()) - 1;
                        int residueNumber2 = Integer.valueOf(tds.get(4).text()) - 1;

                        double plmScore = Double.valueOf(tds.get(6).text());
                        boolean localContact = Math.abs(residueNumber1 - residueNumber2) < 6;

                        if(localContact) {
                            System.out.println("local contact: " + element.text());
                            localPlmScores.get(residueNumber1).add(plmScore);
                            localPlmScores.get(residueNumber2).add(plmScore);
                        } else {
                            System.out.println("long-range contact: " + element.text());
                            longRangePlmScores.get(residueNumber1).add(plmScore);
                            longRangePlmScores.get(residueNumber2).add(plmScore);
                        }
                    });

            return Optional.of(aminoAcids.stream()
                    .map(aminoAcid -> pdbId + ",A," +
                            aminoAcid.getOneLetterCode() + "," +
                            aminoAcid.getResidueIdentifier().getResidueNumber() + "," +
                            (earlyFoldingResidues.contains(aminoAcid) ? "early" : "late") + "," +
                            "local," +
                            StandardFormat.format(localPlmScores.get(aminoAcid.getResidueIndex())
                                    .stream()
                                    .mapToDouble(Double::valueOf)
                                    .average()
                                    .orElse(0.0)) + System.lineSeparator() +
                            pdbId + ",A," +
                            aminoAcid.getOneLetterCode() + "," +
                            aminoAcid.getResidueIdentifier().getResidueNumber() + "," +
                            (earlyFoldingResidues.contains(aminoAcid) ? "early" : "late") + "," +
                            "long-range," +
                            StandardFormat.format(longRangePlmScores.get(aminoAcid.getResidueIndex())
                                    .stream()
                                    .mapToDouble(Double::valueOf)
                                    .average()
                                    .orElse(0.0))
                    )
                    .collect(Collectors.joining(System.lineSeparator())));
        } catch (Exception e) {
            e.printStackTrace();
            return Optional.empty();
        }
    }
}
