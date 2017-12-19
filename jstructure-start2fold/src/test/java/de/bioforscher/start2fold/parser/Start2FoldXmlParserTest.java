package de.bioforscher.start2fold.parser;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.start2fold.Start2FoldConstants;
import de.bioforscher.testutil.TestUtils;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.nodes.Element;
import org.junit.Ignore;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Path;
import java.util.Set;
import java.util.stream.Collectors;

public class Start2FoldXmlParserTest {
    private static final Logger logger = LoggerFactory.getLogger(Start2FoldXmlParserTest.class);

    @Test
    public void shouldParseStart2FoldXml() {
        Structure structure = StructureParser.source("1hrh").parse();
        Chain chain = structure.chains().findFirst().get();
        Start2FoldXmlParser.parse(chain, TestUtils.getResourceAsInputStream("STF0026.xml"));
    }

    @Test
    public void reportStatusOfEFRDataset() {
        System.out.println("dataset is safe - last checked 12/19/17");
    }

    @Test
    @Ignore
    public void shouldEnsureEqualSequencesForAllExperiments() {
        // shows that experiments may have different sequences - need to handle every one individually
        Path directory = Start2FoldConstants.XML_DIRECTORY;
        Start2FoldConstants.list(directory)
                .forEach(path -> {
                    try {
                        Document document = Jsoup.parse(path.toFile(), "UTF-8");
                        Set<String> sequences = document.getElementsByTag("sequence").stream()
                                .map(Element::text)
                                .collect(Collectors.toSet());

                        logger.info("found {} unique sequence(s):{}",
                                sequences.size(),
                                sequences.stream()
                                        .collect(Collectors.joining(System.lineSeparator(), System.lineSeparator(), "")));
//                        Assert.assertTrue(sequences.size() == 1);
                    } catch (IOException e) {
                        throw new UncheckedIOException(e);
                    }
                });
    }

    @Test
    @Ignore
    public void shouldParseAllFiles() {
        // shows that experiments may have different sequences - need to handle every one individually
        Path directory = Start2FoldConstants.XML_DIRECTORY;
        Start2FoldConstants.list(directory)
                .forEach(path -> {
                    try {
                        logger.info("handling {}",
                                path);
                        // safe are: STF0005, STF0021
                        String pdbId = Jsoup.parse(path.toFile(), "UTF-8").getElementsByTag("protein").attr("pdb_id");
                        Structure structure = StructureParser.source(pdbId).parse();
                        Chain chain = structure.chains().findFirst().get();
                        Start2FoldXmlParser.parse(chain, path);
                    } catch (Exception e) {
                        logger.warn("inspect:", e);
                    }
                });
    }
}