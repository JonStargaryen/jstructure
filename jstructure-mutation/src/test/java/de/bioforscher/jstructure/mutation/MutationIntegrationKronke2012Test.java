package de.bioforscher.jstructure.mutation;

import de.bioforscher.jstructure.align.impl.LocalBlastWrapper;
import de.bioforscher.jstructure.feature.evolution.EvolutionaryInformationCalculator;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.mutation.impl.MutationEffectPredictionServiceImpl;
import de.bioforscher.jstructure.mutation.impl.MutationFeatureVectorImpl;
import de.bioforscher.testutil.TestUtils;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Test on the Kronke2012 dataset.
 */
public class MutationIntegrationKronke2012Test {
    private EvolutionaryInformationCalculator evolutionaryInformationCalculator;
    private MutationEffectPredictionService mutationEffectPredictionService;

    @Before
    public void setup() {
        this.evolutionaryInformationCalculator = new EvolutionaryInformationCalculator();
        this.mutationEffectPredictionService = new MutationEffectPredictionServiceImpl();
    }

    @Test
    public void shouldComposeArff() {
        Map<String, List<String[]>> proteinMap = TestUtils.getResourceAsStream("kronke2016/kronke2016.csv")
                .filter(line -> !line.startsWith("original"))
                .map(line -> line.split(","))
                .collect(Collectors.groupingBy(split -> split[18]));

        String output = proteinMap.entrySet()
                .stream()
                .flatMap(entry -> handleLine(entry.getKey(), entry.getValue()))
                .collect(Collectors.joining(System.lineSeparator()));

        System.out.println("@RELATION mutation" + System.lineSeparator() +
                "@ATTRIBUTE pdbId          STRING" + System.lineSeparator() +
                "@ATTRIBUTE mutation       STRING" + System.lineSeparator() +
                "@ATTRIBUTE gutteridge     STRING" + System.lineSeparator() +
                MutationFeatureVectorImpl.toPartialHeader() + System.lineSeparator() +
                "@ATTRIBUTE pssmPenality   NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE shannon        NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE ddg            NUMERIC" + System.lineSeparator() +
                "@ATTRIBUTE class          {effect,tolerated}" + System.lineSeparator() +
                "@DATA");
        System.out.println(output);
    }

    private Stream<String> handleLine(String pdbId, List<String[]> list) {
        try {
            // fetch original structure
            Structure structure = StructureParser.source(pdbId).parse();
            String chainId = structure.select()
                    .residueNumber(Integer.parseInt(list.get(0)[1]))
                    .asAminoAcid()
                    .getParentChain()
                    .getChainIdentifier()
                    .getChainId();

            // isolate chain with mutation - file provides no chainIds - mapping seems safe though
            Structure isolatedStructure = structure.select()
                    .chainName(chainId)
                    .asIsolatedStructure();
            Chain isolatedChain = isolatedStructure.select()
                    .chainName(chainId)
                    .asChain();

            LocalBlastWrapper.PsiBlastResult psiBlastResult = evolutionaryInformationCalculator.composePsiBlastResult(TestUtils.getResourceAsStream("kronke2016/" + pdbId.toLowerCase() + ".out"),
                    TestUtils.getResourceAsStream("kronke2016/" + pdbId.toLowerCase() + ".mtx"));
            evolutionaryInformationCalculator.assignPsiBlastResultToChain(isolatedChain,
                    psiBlastResult);

            MutationJob mutationJob = mutationEffectPredictionService.createMutationJob(pdbId,
                    isolatedChain,
                    psiBlastResult);

            return list.stream()
                    .map(line -> handleLine(mutationJob, line))
                    .filter(Optional::isPresent)
                    .map(Optional::get);
        } catch (Exception e) {
            e.printStackTrace();
            return Stream.empty();
        }
    }

    private Optional<String> handleLine(MutationJob mutationJob, String[] line) {
        try {
            String original = line[0];
            int position = Integer.valueOf(line[1]);
            String mutation = line[2];
            double ddg = Double.valueOf(line[3]);
            double blast = Double.valueOf(line[21]);
            double shannon = Double.valueOf(line[22]);

            AminoAcid originalAminoAcid = mutationJob.getReferenceChain()
                    .select()
                    .residueNumber(position)
                    .asAminoAcid();

            // ensure that selected amino acid is correct
            Assert.assertEquals(original, originalAminoAcid.getOneLetterCode());

            MutationFeatureVector mutationEffectVector = mutationEffectPredictionService.createMutationFeatureVector(mutationJob,
                    originalAminoAcid.getParentChain().getChainIdentifier(),
                    originalAminoAcid.getResidueIdentifier(),
                    AminoAcid.Family.resolveOneLetterCode(mutation));

            String output = mutationEffectVector + "," +
                    blast + "," +
                    shannon + "," +
                    ddg + "," +
                    (ddg < -0.5 ? "effect" : "tolerated");

            System.out.println(output);

            return Optional.of(output);
        } catch (Exception e) {
            e.printStackTrace();
            return Optional.empty();
        }
    }
}
