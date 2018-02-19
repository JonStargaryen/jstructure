package de.bioforscher.jstructure.feature.uniprot.homologous;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.identifier.ChainIdentifier;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.kraken.interfaces.uniprot.features.FeatureSequence;
import uk.ac.ebi.kraken.interfaces.uniprot.features.FeatureType;
import uk.ac.ebi.kraken.interfaces.uniprot.features.MutagenFeature;

import java.util.List;
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
        Structure protein = StructureParser.fromPdbId("2ts1").parse();
        uniProtHomologyAnnotator.uniProtBlastService.start();
        uniProtHomologyAnnotator.process(protein);
        uniProtHomologyAnnotator.uniProtBlastService.stop();

        List<MutagenFeature> mutagenFeatures = protein.chainsWithAminoAcids()
                .flatMap(Chain::aminoAcids)
                .flatMap(aminoAcid -> aminoAcid.getFeature(UniProtFeatureContainer.class)
                        .getFeatures(FeatureType.MUTAGEN)
                        .stream())
                .filter(MutagenFeature.class::isInstance)
                .map(MutagenFeature.class::cast)
                .collect(Collectors.toList());
        Assert.assertFalse("did not find any mutagen features", mutagenFeatures.isEmpty());
        mutagenFeatures.forEach(mutagenFeature -> {
            System.out.println(mutagenFeature.getMutagenReport().getValue());
            System.out.println(mutagenFeature.getOriginalSequence().getValue() + " -> " +
                    mutagenFeature.getAlternativeSequences()
                            .stream()
                            .map(FeatureSequence::getValue)
                            .collect(Collectors.toList()));
        });

        System.out.println("homologous protein chains:");
        List<ChainIdentifier> chains = protein.chainsWithAminoAcids()
                .flatMap(chain -> chain.getFeature(UniProtHomologousEntryContainer.class)
                        .getHomologousChains()
                        .stream())
                .collect(Collectors.toList());
        Assert.assertTrue("did not find homologous chains", chains.size() > 1);
        chains.forEach(System.out::println);
    }
}