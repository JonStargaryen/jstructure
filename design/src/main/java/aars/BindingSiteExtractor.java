package aars;

import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.family.AminoAcidFamily;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import de.bioforscher.jstructure.parser.ProteinParser;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Extract binding sites based on Chris' data.
 * Created by bittrich on 2/8/17.
 */
public class BindingSiteExtractor {
    public static void main(String[] args) {
        new BindingSiteExtractor().run();
    }

    private static final List<String> classNames = Stream.of("C1", "C2").collect(Collectors.toList());
    private static final List<String> modeNames = Stream.of("1", "2_aa_part", "2_axp_part", "3", "4").collect(Collectors.toList());

    private void run() {
        for(String className : classNames) {
            for(String modeName : modeNames) {
                Path path = Paths.get(AARSConstants.BINDING_SITE_PATH + className + "/ligand_based/per_type/contacts_mode_" + modeName);
                AARSConstants.lines(path).forEach(line -> handleLine(line, className, modeName));
            }
        }
    }

    private void handleLine(String line, String className, String modeName) {
        System.out.println(line);

        // map to amino acid family
        String aminoAcid = Stream.of(AminoAcidFamily.values())
                .filter(aminoAcidFamily -> aminoAcidFamily.name().equalsIgnoreCase(line.split(":")[0]))
                .findFirst()
                .map(AminoAcidFamily::getThreeLetterCode)
                .orElse("Pyr");

        int[] residueNumbers = Pattern.compile(", ").splitAsStream(line.split("\\[")[1].split("]")[0])
                .mapToInt(Integer::valueOf)
                .toArray();

        AARSConstants.lines(Paths.get(AARSConstants.MAIN_TABLE_CATALYTIC_PATH))
                // describing this amino acid
                .filter(mainTableLine -> mainTableLine.split(",")[1].equalsIgnoreCase(aminoAcid))
                // representative
                .filter(mainTableLine -> mainTableLine.endsWith("1"))
                .forEach(mainTableLine -> {
                    String[] split = mainTableLine.split(",");
                    String pdbId = split[2];
                    String chainId = split[13];
                    Protein protein = ProteinParser.source(AARSConstants.RENUMBERED_STRUCTURES_PATH + className + "/renumbered_structures/" + pdbId + "_renum.pdb").parse();

                    String output = "HEADER    LIGASE/RNA                              15-JUL-99   " + pdbId.toUpperCase() + "              " + System.lineSeparator() + Selection.on(protein)
                            .chainName(chainId)
                            .aminoAcids()
                            .residueNumber(residueNumbers)
                            .asGroupContainer()
                            .composePDBRecord();
                    AARSConstants.write(Paths.get(AARSConstants.BINDING_SITE_PATH + className + "/ligand_based/per_type/"
                            + modeName +  "/" + aminoAcid.substring(0, 1).toUpperCase() + aminoAcid.substring(1, 3).toLowerCase() + "/" +
                            pdbId + "_" + chainId + "_renum.pdb"), output.getBytes());
                });
    }
}
