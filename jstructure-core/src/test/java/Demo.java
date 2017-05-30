import de.bioforscher.jstructure.feature.loopfraction.LoopFraction;
import de.bioforscher.jstructure.model.feature.FeatureProviderRegistry;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.aminoacid.Alanine;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Test;

/**
 * Some code fragments demonstrating some functions of the library.
 * Created by S on 03.11.2016.
 */
public class Demo {
    @Test
    public void demo() {
        // fetch protein from the PDB
        Protein protein = ProteinParser.source("1brr").parse();

        // compute the fraction of alanine in the protein
        double alanineRatio = protein.select()
                .groupName(Alanine.THREE_LETTER_CODE)
                .asFilteredGroups()
                .count() / (double) protein.getSize() * 100.0;
        System.out.printf("alanine ratio: %3.2f%%", alanineRatio);

        // compute some features - e.g. compute the loop fraction
        FeatureProviderRegistry.resolve(LoopFraction.class).process(protein);
    }
}
