package de.bioforscher.jstructure.feature.uniprot.homologous;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.ProteinParser;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.kraken.interfaces.uniprot.features.FeatureSequence;
import uk.ac.ebi.kraken.interfaces.uniprot.features.MutagenFeature;

import java.util.stream.Collectors;

/**
 * Test the UniProt annotator based on homology.
 * Created by bittrich on 7/10/17.
 */
public class UniProtHomologyAnnotatorTest {
    private UniProtHomologyAnnotator uniProtHomologyAnnotator;

    @Before
    public void setup() {
        uniProtHomologyAnnotator = new UniProtHomologyAnnotator();
    }

    @Test
    public void shouldInferInformationFromHomologousProteins() {
        Protein protein = ProteinParser.source("2ts1").parse();
        uniProtHomologyAnnotator.uniProtBlastService.start();
        uniProtHomologyAnnotator.process(protein);
        uniProtHomologyAnnotator.uniProtBlastService.stop();

        protein.chainsWithAminoAcids()
                .flatMap(Chain::aminoAcids)
                .filter(aminoAcid -> aminoAcid.getFeatureContainer()
                        .getFeature(UniProtFeatureContainer.class)
                        .getFeatures()
                        .values()
                        .stream()
                        .anyMatch(feature -> feature instanceof MutagenFeature))
                .forEach(aminoAcid -> {
                    System.out.println(aminoAcid);
                    aminoAcid.getFeatureContainer()
                            .getFeature(UniProtFeatureContainer.class)
                            .getFeatures()
                            .entrySet()
                            .stream()
                            .filter(entry -> entry.getValue() instanceof MutagenFeature)
                            .forEach(entry -> {
                                System.out.println("source: " + entry.getKey());
                                MutagenFeature mutagenFeature = (MutagenFeature) entry.getValue();
                                System.out.println(mutagenFeature.getMutagenReport().getValue());
                                System.out.println(mutagenFeature.getOriginalSequence().getValue() + " -> " +
                                        mutagenFeature.getAlternativeSequences()
                                                .stream()
                                                .map(FeatureSequence::getValue)
                                                .collect(Collectors.toList()));
                            });
                });
    }
}