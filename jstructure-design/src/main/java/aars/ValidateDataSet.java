package aars;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

/**
 * Validate the main table by the following criteria:
 * <ul>
 *     <li>resolution is correct</li>
 *     <li>organism mapping is correct</li>
 *     <li>wild-type annotation</li>
 *     <li>ligands can be found at specified position and map to the correct name</li>
 * </ul>
 * Created by S on 12.05.2017.
 */
public class ValidateDataSet {
    private static final Path mainTablePath = Paths.get(System.getProperty("user.home"),
            "/git/aars_data/main_table/main_table_curated_v01.csv");

    public static void main(String[] args) throws IOException {
        List<String> errors = Files.lines(mainTablePath)
                // skip header line
                .filter(line -> !line.startsWith("identifier"))
                .map(ValidateDataSet::handleLine)
                .filter(line -> !line.isEmpty())
                .collect(Collectors.toList());

        System.err.println("errors:");
        errors.forEach(System.err::println);
    }

    private static String handleLine(String line) {
        try {
            handleLineChecked(line);
            return "";
        } catch(Exception e) {
            System.err.println(e.getLocalizedMessage());
            return line + ";" + e.getLocalizedMessage();
        }
    }

    private static void handleLineChecked(String line) {
        System.out.println(line);
        String[] split = line.split(";");
        String pdbId = split[0].split("_")[0];
        String chainId = split[0].split("_")[1];

        Protein protein = ProteinParser.source(pdbId)
                .cachedMode()
                .parse();

        // ensure chain is present
        Chain chain = protein
                .select()
                .chainName(chainId)
                .asChain();

        String ligandName = split[7];

        if(!ligandName.equals("x")) {
            // check if ligand is present at defined position
            Optional<Group> ligand = protein.select()
                    .chainName(split[8])
                    .residueNumber(Integer.valueOf(split[9]))
                    .asOptionalGroup();

            if(!ligand.isPresent() || !ligandName.equals(ligand.get().getThreeLetterCode())) {
                throw new IllegalArgumentException("did not find ligand: " + ligandName + " at " + split[8] + "-" +
                        split[9]);
            }

            // check that amino acid ligands map to the correct binary mode
            if(ligand.get().isAminoAcid()) {
                if(!split[10].equals("1") ||!split[11].equals("0") || !split[12].equals("1")) {
                    throw new IllegalArgumentException("binary mode string does not exception for plain amino acid, found: [" + split[10] + ", " + split[11] + ", " + split[12] + "]");
                }
            }
        }
    }
}