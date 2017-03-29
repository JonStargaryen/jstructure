package aars;

import de.bioforscher.jstructure.model.structure.family.GroupInformation;
import de.bioforscher.jstructure.parser.CIFParser;
import de.bioforscher.jstructure.parser.sifts.SiftsParser;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.nodes.Element;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.DoubleSummaryStatistics;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Enrich contour plot with additional data of the proteins function (maybe UniProt?) and ligand type (long name from
 * CIF-file).
 * Created by bittrich on 3/29/17.
 */
public class ContourPlotEnricher {
    public static void main(String[] args) throws IOException {
//        deriveStatistics();

        extractRegionOfPlot().forEach(ContourPlotEnricher::determineMapping);
        Files.write(Paths.get("/home/bittrich/git/aars_analysis/data/geometry/saltbridges_function.tsv"), output.toString().getBytes());
    }

    private static final StringBuilder output = new StringBuilder();

    /**
     * Map a given container to protein and ligand function.
     * @param ofInterest the container full of gold
     */
    private static void determineMapping(OfInterest ofInterest) {
        SiftsParser siftsParser = new SiftsParser();
        List<String> uniProtIds = siftsParser.mapToUniProt(ofInterest.pdbId, ofInterest.chainId);

        System.out.println(ofInterest.pdbId + "_" + ofInterest.chainId + " : " + ofInterest.ligand.getName());
        uniProtIds.forEach(uniProtId -> {
            try {
                Document document = Jsoup.connect("http://www.uniprot.org/uniprot/" + uniProtId + ".xml").get();
                List<String> keywords = document.getElementsByTag("keyword").stream()
                        .map(Element::text)
                        .filter(ContourPlotEnricher::filterKeywords)
                        .collect(Collectors.toList());
                System.out.println(keywords);
                output.append(ofInterest).append("\t").append(ofInterest.ligand.getName().replace("\"", "")).append("\t").append(uniProtId).append("\t").append(keywords).append(System.lineSeparator());
            } catch (IOException e) {
                System.err.println("failed to load UniProt entry '" + uniProtId + "'");
            }
        });
    }

    private static final List<String> forbiddenKeywords = Stream.of("3D-structure",
            "Alternative splicing",
            "Complete proteome",
            "Cytoplasm",
            "Direct protein sequencing",
            "Polymorphism",
            "Reference proteome").collect(Collectors.toList());

    private static boolean filterKeywords(String keyword) {
        return !forbiddenKeywords.contains(keyword);
    }

    /**
     * Gather all occurrences of salt-bridge interaction similar to aaRS tweezers.
     * @return all relevant containers
     */
    private static Stream<OfInterest> extractRegionOfPlot() {
        output.append("pairId\tcaDistance\tlastCarbonDistance\tbackboneAngle\tsideChainAngle\tligand\tuniProtId\tkeywords").append(System.lineSeparator());
        return AARSConstants.lines(Paths.get("/home/bittrich/git/aars_analysis/data/geometry/saltbridges_geometry_phosphate.tsv"))
                .filter(line -> !line.startsWith("pairId"))
                .map(line -> line.split("\t"))
                .filter(ContourPlotEnricher::isInRegionOfInterest)
//                .limit(10)
                .map(OfInterest::new)
                .filter(OfInterest::isSane);
    }

    /**
     * True iff similar to observed arg tweezers conformations.
     * @param split data to decide on
     * @return true when relevant
     */
    private static boolean isInRegionOfInterest(String[] split) {
        double minDistance = 13.39;
        double maxDistance = 16.03;
        double minAngle = 54.15;
        double maxAngle = 111.08;

        try {
            double distance = Double.valueOf(split[1]);
            double angle = Double.valueOf(split[4]);
            return distance > minDistance && distance < maxDistance && angle > minAngle && angle < maxAngle;
        } catch (ArrayIndexOutOfBoundsException e) {
            // some lines to not contain information
            return false;
        }

    }

    static class OfInterest {
        String pdbId, chainId, ligandId;
        GroupInformation ligand;
        boolean sane = true;
        String originalLine;

        OfInterest(String[] split) {
            originalLine = Stream.of(split).collect(Collectors.joining("\t"));
            String[] idSplit = split[0].split(":");
            String[] rangeSplit = idSplit[3].split("-");
            this.pdbId = idSplit[0];
            this.chainId = idSplit[2];
            this.ligandId = idSplit[1];
            this.ligand = CIFParser.parseLigandInformation(ligandId);
            if(rangeSplit[1].charAt(rangeSplit[1].length() - 1) != rangeSplit[2].charAt(rangeSplit[2].length() - 1)) {
                System.err.println("not matching chain ids");
                this.sane = false;
            }
        }

        boolean isSane() {
            return sane;
        }

        @Override
        public String toString() {
            return originalLine;
        }
    }

    /**
     * Gather statistics on which observations of salt-bridges are similar to arg tweezers.
     */
    private static void deriveStatistics() {
        List<String[]> lines = AARSConstants.lines(Paths.get("/home/bittrich/git/aars_analysis/data/geometry/argtweezer_geometry.tsv"))
                .filter(line -> !line.startsWith("id"))
                .map(line -> line.split("\t"))
                .filter(split -> split[3].equals("A"))
                .collect(Collectors.toList());

        DoubleSummaryStatistics distanceStatistics = lines.stream()
                .map(line -> line[1])
                .mapToDouble(Double::valueOf)
                .summaryStatistics();
        DoubleSummaryStatistics angleStatistics = lines.stream()
                .map(line -> line[2])
                .mapToDouble(Double::valueOf)
                .summaryStatistics();

        System.out.println("min-dist: " + distanceStatistics.getMin());
        System.out.println("max-dist: " + distanceStatistics.getMax());
        System.out.println("min-angle: " + angleStatistics.getMin());
        System.out.println("max-angle: " + angleStatistics.getMax());
    }
}
