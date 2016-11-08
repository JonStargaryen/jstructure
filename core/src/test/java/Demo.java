import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.Residue;
import de.bioforscher.jstructure.model.structure.filter.AminoAcidFilter;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Test;

import java.util.Collections;
import java.util.function.Predicate;

import static de.bioforscher.jstructure.model.structure.AminoAcid.ALANINE;

/**
 * Some code fragments demonstrating some functions of the library.
 * Created by S on 03.11.2016.
 */
public class Demo {
    enum FeatureNames {
        ALANINE_RATIO
    }

    @Test
    public void shouldPrintAllResiduesInContact() {
        // fetch/parse structure by id
        Protein protein = ProteinParser.parseProteinById("1brr");

        // print all getResidue pairs whose distance is less than 8.0 A
        protein.residuePairsInContact(8.0)
               .forEach(System.out::println);

        // custom filters
        Predicate<Residue> alanineFilter = new AminoAcidFilter(Collections.singletonList(ALANINE));
        // print coordinates of all alanines
        protein.residues()
               .filter(alanineFilter)
               .map(Residue::composePDBRecord)
               .forEach(System.out::println);

        // count all alanines
        double alanineRatio = protein.residues()
                .filter(alanineFilter)
                .count() / (double) protein.getSize() * 100.0;

        // store count
        protein.setFeature(FeatureNames.ALANINE_RATIO, alanineRatio);

        // retrieve it
        System.out.printf("alanine ratio: %3.2f%%", protein.getDoubleFeature(FeatureNames.ALANINE_RATIO));
    }
}
