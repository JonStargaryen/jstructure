package de.bioforscher.jstructure.si.explorer;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.efr.model.si.ContactStructuralInformation;
import de.bioforscher.jstructure.efr.model.si.ResidueStructuralInformation;
import de.bioforscher.jstructure.efr.model.SecondaryStructureElement;
import de.bioforscher.jstructure.efr.model.Start2FoldResidueAnnotation;
import de.bioforscher.jstructure.efr.parser.EvolutionaryCouplingParser;
import de.bioforscher.jstructure.efr.parser.Start2FoldXmlParser;
import de.bioforscher.jstructure.efr.parser.StructuralInformationParserService;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.testutil.TestUtils;
import org.jsoup.Jsoup;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.List;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class ExplorerChain {
    private static final Logger logger = LoggerFactory.getLogger(ExplorerChain.class);
    private static final StructuralInformationParserService STRUCTURAL_INFORMATION_PARSER_SERVICE =
            StructuralInformationParserService.getInstance();
    private static final Pattern PATTERN = Pattern.compile(",");

    private final String stfId;
    private final String pdbId;
    private final List<Integer> experimentIds;
    private final String uniProtId;
    private final List<Integer> earlyResidueNumbers;
    private final List<Integer> functionalResidueNumbers;
    private final String pdbRepresentation;
    private final List<ContactStructuralInformation> contacts;
    private final List<ResidueStructuralInformation> residues;
    private final String title;
    private final String sequence;
    private final List<SecondaryStructureElement> secondaryStructureElements;
    private final String averageRmsd;
    private final String averageTmScore;
    private final String averageQ;

    public ExplorerChain(String[] lineSplit) {
        // parse information of list
        this.stfId = lineSplit[0];
        this.pdbId = lineSplit[1];
        this.experimentIds = parseIntegerList(lineSplit[2]);
        this.uniProtId = lineSplit[4];
        this.functionalResidueNumbers = parseIntegerList(lineSplit[5]);

        logger.info("creating chain {} - {}",
                this.stfId,
                this.pdbId);

        // parse chain & compute features
        Chain chain = StructureParser.fromInputStream(TestUtils.getResourceAsInputStream("data/pdb/" + stfId + ".pdb"))
                .parse()
                .getFirstChain();
        this.pdbRepresentation = chain.getPdbRepresentation();

        Start2FoldXmlParser.parseSpecificExperiment(chain,
                TestUtils.getResourceAsInputStream("data/xml/" + stfId + ".xml"),
                experimentIds);

        List<AminoAcid> aminoAcids = chain.aminoAcids()
                .collect(Collectors.toList());
        List<AminoAcid> earlyFoldingResidues = aminoAcids.stream()
                .filter(aminoAcid -> aminoAcid.getFeature(Start2FoldResidueAnnotation.class).isEarly())
                .collect(Collectors.toList());
        this.earlyResidueNumbers = earlyFoldingResidues.stream()
                .map(AminoAcid::getResidueIdentifier)
                .map(ResidueIdentifier::getResidueNumber)
                .collect(Collectors.toList());

        EvolutionaryCouplingParser.parseHotSpotFile(chain,
                TestUtils.getResourceAsInputStream("data/coupling/" + stfId + "_hs.html"));
        this.contacts = STRUCTURAL_INFORMATION_PARSER_SERVICE.parseContactStructuralInformationFile(TestUtils.getResourceAsInputStream("data/raw/" + stfId + ".out"),
                earlyFoldingResidues);
        this.residues = STRUCTURAL_INFORMATION_PARSER_SERVICE.composeResidueStructuralInformation(aminoAcids,
                earlyFoldingResidues,
                contacts);
        try {
            EvolutionaryCouplingParser.parsePlmScore(contacts,
                    Jsoup.parse(TestUtils.getResourceAsInputStream("data/coupling/" + stfId + "_ec.html"), "UTF-8", ""));
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
        double[] rawBaselinePerformance = STRUCTURAL_INFORMATION_PARSER_SERVICE.parseAveragePerformance(TestUtils.getResourceAsInputStream("data/raw/" + stfId + ".out"));
        this.averageRmsd = StandardFormat.format(rawBaselinePerformance[0]);
        this.averageTmScore = StandardFormat.format(rawBaselinePerformance[1]);
        this.averageQ = StandardFormat.format(rawBaselinePerformance[2]);

        this.title = chain.getParentStructure().getTitle();
        this.sequence = chain.getAminoAcidSequence();

        this.secondaryStructureElements = SecondaryStructureElement.of(chain);
    }

    private static List<Integer> parseIntegerList(String element) {
        return PATTERN.splitAsStream(element.replace("[", "").replace("]", ""))
                .flatMap(ExplorerChain::parseInteger)
                .collect(Collectors.toList());
    }

    private static Stream<Integer> parseInteger(String value) {
        if(value.contains("-")) {
            String[] split = value.split("-");
            return IntStream.range(Integer.valueOf(split[0]), Integer.valueOf(split[1]))
                    .boxed();
        } else {
            return Stream.of(Integer.valueOf(value));
        }
    }

    public String getStfId() {
        return stfId;
    }

    public String getPdbId() {
        return pdbId;
    }

    public List<Integer> getExperimentIds() {
        return experimentIds;
    }

    public String getUniProtId() {
        return uniProtId;
    }

    public List<Integer> getEarlyResidueNumbers() {
        return earlyResidueNumbers;
    }

    public List<Integer> getFunctionalResidueNumbers() {
        return functionalResidueNumbers;
    }

    public String getPdbRepresentation() {
        return pdbRepresentation;
    }

    public List<ContactStructuralInformation> getContacts() {
        return contacts;
    }

    public List<ResidueStructuralInformation> getResidues() {
        return residues;
    }

    public String getAverageRmsd() {
        return averageRmsd;
    }

    public String getAverageTmScore() {
        return averageTmScore;
    }

    public String getAverageQ() {
        return averageQ;
    }

    public String getTitle() {
        return title;
    }

    public String getSequence() {
        return sequence;
    }

    public List<SecondaryStructureElement> getSecondaryStructureElements() {
        return secondaryStructureElements;
    }
}
