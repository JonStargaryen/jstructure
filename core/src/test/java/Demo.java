import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.family.AminoAcidFamily;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Test;

/**
 * Some code fragments demonstrating some functions of the library.
 * Created by S on 03.11.2016.
 */
public class Demo {
    enum FeatureNames {
        ALANINE_RATIO
    }

    @Test
    public void demo() {
        // fetch/parse structure by id
        Protein protein = ProteinParser.parseProteinById("1brr");

        // print all getResidue pairs whose distance is less than 8.0 A
        //TODO pairsOn(Container) seems to be horribly slow - move to for-loops or drop all together?
        Selection.pairsOn(protein)
                .alphaCarbonDistance(8.0)
                .asFilteredGroupPairs()
                .forEach(System.out::println);

        // print coordinates of all alanines
        Selection.on(protein)
                .aminoAcids(AminoAcidFamily.ALANINE)
                .asFilteredGroups()
                .map(Group::composePDBRecord)
                .forEach(System.out::println);

        // count all alanines
        double alanineRatio = Selection.on(protein)
                .aminoAcids(AminoAcidFamily.ALANINE)
                .asFilteredGroups()
                .count() / (double) protein.getSize() * 100.0;

        // store count
        protein.setFeature(FeatureNames.ALANINE_RATIO, alanineRatio);

        // retrieve it
        System.out.printf("alanine ratio: %3.2f%%", protein.getFeatureAsDouble(FeatureNames.ALANINE_RATIO));
    }
}
