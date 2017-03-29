package aars;

import de.bioforscher.jstructure.parser.sifts.SiftsParser;
import de.bioforscher.jstructure.parser.uniprot.UniProtAnnotationContainer;
import de.bioforscher.jstructure.parser.uniprot.UniProtAnnotator;
import de.bioforscher.jstructure.parser.uniprot.UniProtMutagenesisSite;
import de.bioforscher.jstructure.parser.uniprot.UniProtNaturalVariant;

import java.nio.file.Paths;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Assign uniprot id to all entries in the main table.
 * Created by bittrich on 3/28/17.
 */
public class UniProtMutationAnnotator {
    private static final SiftsParser siftsParser = new SiftsParser();
    private static final UniProtAnnotator uniProtAnnotator = new UniProtAnnotator();
    private static final String DEL = "\t";

    public static void main(String[] args) {
        System.out.println("pdbId" + DEL + "chainId" + DEL + "uniprotId" + DEL + "type" + DEL + "pos" + DEL + "original" + DEL + "variant" + DEL + "description");
        AARSConstants.lines(Paths.get(AARSConstants.MAIN_TABLE_CATALYTIC_PATH))
                .filter(line -> !line.startsWith("class"))
                .map(line -> line.split(","))
                .forEach(UniProtMutationAnnotator::process);
    }

    private static void process(String[] split) {
        String pdbId = split[2];
        String chainId = split[13];
        List<String> uniProtId = siftsParser.mapToUniProt(pdbId, chainId);
//        System.out.println(pdbId + "_" + chainId + " -> " + uniProtId);

        uniProtId.stream()
                .map(uniProtAnnotator::process)
                .map(container -> renderContainer(pdbId, chainId, container))
                .filter(string -> !string.isEmpty())
                .forEach(System.out::println);
    }

    private static String renderContainer(String pdbId, String chainId, UniProtAnnotationContainer container) {
        List<String> mutations = container.getMutagenesisSites().stream()
                .map(UniProtMutationAnnotator::renderMutation)
                .collect(Collectors.toList());

        List<String> variants = container.getNaturalVariants().stream()
                .map(UniProtMutationAnnotator::renderVariation)
                .collect(Collectors.toList());

        return Stream.of(mutations, variants)
                .flatMap(Collection::stream)
                .map(string -> pdbId + DEL + chainId + DEL + container.getUniProtId() + DEL + string)
                .collect(Collectors.joining(System.lineSeparator()));
    }

    private static String renderVariation(UniProtNaturalVariant naturalVariant) {
        return "V" + DEL +
                naturalVariant.getPosition() + DEL +
                naturalVariant.getOriginal() + DEL +
                naturalVariant.getVariation() + DEL +
                naturalVariant.getDescription();
    }

    private static String renderMutation(UniProtMutagenesisSite mutagenesisSite) {
        return "M" + DEL +
                mutagenesisSite.getPosition() + DEL +
                mutagenesisSite.getOriginal() + DEL +
                mutagenesisSite.getVariation() + DEL +
                mutagenesisSite.getDescription();
    }
}
