package de.bioforscher.jstructure.mutation.old.impl;

import de.bioforscher.jstructure.feature.uniprot.UniProtBridge;
import de.bioforscher.jstructure.model.feature.FeatureProviderRegistry;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.identifier.ChainIdentifier;
import de.bioforscher.jstructure.mutation.old.MutationDescriptor;
import de.bioforscher.jstructure.mutation.old.MutationJob;
import de.bioforscher.testutil.TestUtils;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import uk.ac.ebi.kraken.interfaces.uniprot.PrimaryUniProtAccession;
import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry;
import uk.ac.ebi.uniprot.dataservice.client.uniprot.UniProtQueryBuilder;
import uk.ac.ebi.uniprot.dataservice.client.uniprot.UniProtService;

import java.io.FileWriter;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Test for the function of the mutation effect predictor.
 * Created by bittrich on 7/11/17.
 */
public class MutationEffectPredictionServiceImplTest {
    private static final Logger logger = LoggerFactory.getLogger(MutationEffectPredictionServiceImplTest.class);
    private MutationEffectPredictionServiceImpl mutationEffectPredictor;
    private UniProtService uniProtService;

    @Before
    public void setup() {
        mutationEffectPredictor = new MutationEffectPredictionServiceImpl();
        uniProtService = UniProtBridge.getInstance().getUniProtService();
        // access the feature registry once, so logs are 'in-order'
        FeatureProviderRegistry.getRegisteredFeatureProviders();
    }

    @Test
    public void shouldBlastForHomologousSequenceAndStructures() {
        MutationJob mutationJob = new MutationJobImpl("CYC32_DESDN",
                "ETFEIPESVTMSPKQFEGYTPKKGDVTFNHASHMDIACQQCHHTVPDTYTIESCMTEGCHDNIKERTEISSVYRTFHTTKDSEKSCVGCHRELKRQGPSDAPLACNSCHVQ");
        mutationEffectPredictor.createMultiSequenceAlignment(mutationJob);

        System.out.println("homologous sequences:");
        System.out.println(mutationJob.getHomologousEntryContainer()
                .getUniProtEntries()
                .stream()
                .map(UniProtEntry::getPrimaryUniProtAccession)
                .map(PrimaryUniProtAccession::getValue)
                .collect(Collectors.toList()));
        Assert.assertTrue(mutationJob.getHomologousEntryContainer().getUniProtEntries().size() >= 8);

        System.out.println("homologous pdb chains:");
        System.out.println(mutationJob.getHomologousPdbChains()
                .stream()
                .map(Chain::getChainIdentifier)
                .map(ChainIdentifier::getFullName)
                .collect(Collectors.toList()));
        Assert.assertTrue(mutationJob.getHomologousPdbChains().size() > 0);
    }

    @Test
    public void shouldComposeMutationDescriptor() {
        MutationJob mutationJob = new MutationJobImpl("CYC32_DESDN",
                "ETFEIPESVTMSPKQFEGYTPKKGDVTFNHASHMDIACQQCHHTVPDTYTIESCMTEGCHDNIKERTEISSVYRTFHTTKDSEKSCVGCHRELKRQGPSDAPLACNSCHVQ");
        mutationEffectPredictor.createMultiSequenceAlignment(mutationJob);
        System.out.println(mutationJob.composeMutationDescriptor(73, AminoAcid.Family.GLUTAMIC_ACID));
    }

    @Test
    @Ignore
    public void shouldComposeSchaefer2012SequenceDataset() throws IOException {
        // key: id - value: full line
        Map<String, List<String[]>> sequenceMap = TestUtils.getResourceAsStream("schaefer2012/mutants_function.txt")
                .filter(line -> !line.startsWith("ID"))
                .map(line -> line.split(","))
                .limit(50)
                .collect(Collectors.groupingBy(split -> split[0]));
        FileWriter fileWriter = new FileWriter("/home/bittrich/schaefer2012-sequence.arff", true);

        sequenceMap.keySet()
                .stream()
                .flatMap(id -> handleSequenceBin(sequenceMap, id))
                .forEach(line -> {
                    try {
                        System.out.println(line);
                        fileWriter.write(line + System.lineSeparator());
                    } catch (IOException e) {
                        throw new UncheckedIOException(e);
                    }
                });
    }

    @Test
    @Ignore
    public void shouldComposeSchaefer2012StructureDataset() throws IOException {
        // key: chainId - value: full line
        Map<String, List<String[]>> chainMap = TestUtils.getResourceAsStream("schaefer2012/mutants_stability.txt")
                .filter(line -> !line.startsWith("ID"))
                .map(line -> line.split(","))
                .collect(Collectors.groupingBy(split -> split[0]));

        Path outputPath = Paths.get("/home/bittrich/schaefer2012-structure.arff");
        boolean writeHeader = !Files.exists(outputPath);
        FileWriter fileWriter = new FileWriter(outputPath.toFile(), true);
        if(writeHeader) {
            fileWriter.append(getHeaderSchaefer2012());
        }

        chainMap.keySet()
                .stream()
                .flatMap(id -> handleStructureBin(chainMap, id))
                .forEach(line -> {
                    try {
                        System.out.println(line);
                        fileWriter.write(line + System.lineSeparator());
                        fileWriter.flush();
                    } catch (IOException e) {
                        throw new UncheckedIOException(e);
                    }
                });

        fileWriter.close();
    }

    private Stream<String> handleSequenceBin(Map<String, List<String[]>> sequenceMap, String id) {
        try {
            // old sequence identifiers can be retrieved via http://www.uniprot.org/uniprot/%s.fasta
            String sequence = uniProtService.getEntries(UniProtQueryBuilder.id(id)).getFirstResult().getSequence().getValue();

            MutationJob mutationJob = new MutationJobImpl(id, sequence);
            mutationEffectPredictor.createMultiSequenceAlignment(mutationJob);

            return sequenceMap.get(id)
                    .stream()
                    .map(line -> handleSequenceLine(mutationJob, line))
                    .filter(Optional::isPresent)
                    .map(Optional::get);
        } catch (Exception e) {
            // will fail to resolve some ids
            logger.warn("failed to resolve sequence for {}", id);
            return Stream.empty();
        }
    }

    @Test
    public void moreErrorCases() {
        System.out.println("1oox, 1uo7 pdb missing");
        System.out.println("1bni_c - position 110 missing - java.lang.IndexOutOfBoundsException: Index: 109, Size: 108 - composeMutationDescriptor(MutationJobImpl.java:178)");
    }

    private Stream<String> handleStructureBin(Map<String, List<String[]>> chainMap, String id) {
        try {
            Chain chain = StructureParser.source(id.substring(0, 4))
                    .minimalParsing(true)
                    .parse()
                    .select()
                    .chainName(id.substring(4))
                    .asChain();
            String sequence = chain.getAminoAcidSequence();

            MutationJob mutationJob = new MutationJobImpl(id, sequence, chain);
            mutationEffectPredictor.createMultiSequenceAlignment(mutationJob);

            return chainMap.get(id)
                    .stream()
                    .map(line -> handleStructureLine(mutationJob, line))
                    .filter(Optional::isPresent)
                    .map(Optional::get);
        } catch (Exception e) {
            // will fail to resolve some ids
            logger.warn("computation failed for {}", id, e);
            return Stream.empty();
        }
    }

    private Optional<String> handleSequenceLine(MutationJob mutationJob, String[] line) {
        try {
            int position = Integer.valueOf(line[1]);
            AminoAcid.Family targetAminoAcid = AminoAcid.Family.resolveOneLetterCode(line[2]);
            String hasEffect = line[3].equals("1") ? "effect" : "tolerated";
            MutationDescriptor mutationDescriptor = mutationJob.composeMutationDescriptor(position, targetAminoAcid);
            return Optional.of(mutationDescriptor.asArffString() + "," + hasEffect);
        } catch (Exception e) {
            logger.warn("failed to handle line {}", Arrays.toString(line), e);
            return Optional.empty();
        }
    }

    private Optional<String> handleStructureLine(MutationJob mutationJob, String[] line) {
        try {
            int position = Integer.valueOf(line[1]);
            AminoAcid.Family targetAminoAcid = AminoAcid.Family.resolveOneLetterCode(line[2]);
            String hasEffect = Math.abs(Double.valueOf(line[3])) > 1.0 ? "effect" : "tolerated";
            MutationDescriptor mutationDescriptor = mutationJob.composeMutationDescriptor(position, targetAminoAcid);
            return Optional.of(mutationDescriptor.asArffString() + "," + hasEffect);
        } catch (Exception e) {
            logger.warn("failed to handle line {}", Arrays.toString(line), e);
            return Optional.empty();
        }
    }

    @Test
    @Ignore
    public void printHeaderSchaefer2012() {
        System.out.println(getHeaderSchaefer2012());
    }

    private String getHeaderSchaefer2012() {
        return "@RELATION mutation" + System.lineSeparator() +
                "@ATTRIBUTE chainId        STRING" + System.lineSeparator() +
                "@ATTRIBUTE mutation       STRING" + System.lineSeparator() +
                "@ATTRIBUTE oldFunctional  NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE newFunctional  NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE gutterChanged  NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE aaMaxRasaDelta NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE mutationFrq    NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE mutationObs    NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE deletionFrq    NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE deletionObs    NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE mutToSimGutFrq NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE mutToSimGutObs NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE mutToTargetFrq NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE mutToTargetObs NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE bindSiteFrq    NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE bindSiteObs    NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE interesSiteFrq NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE interesSiteObs NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE mutagenFrq     NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE mutagenObs     NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE variantFrq     NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE variantObs     NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE disulfFrq      NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE disulfObs      NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE disulfCapLost  NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE envRasaDelta   NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE envEnerDelta   NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE envLigConDelta NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE envLoopDelta   NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE aaRasaDelta    NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE aaEnerDelta    NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE aaLigConDelta  NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE aaLoopDelta    NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE class          {effect,tolerated}" + System.lineSeparator() +
                "@DATA" + System.lineSeparator();
    }
}