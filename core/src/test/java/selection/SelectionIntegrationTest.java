package selection;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Before;
import org.junit.Test;

/**
 * Validates the behaviour of the selection API.
 * Created by S on 21.11.2016.
 */
public class SelectionIntegrationTest {
    private Protein protein;

    @Before
    public void setup() {
        protein = ProteinParser.parseProteinById("1brr");
    }

    @Test
    public void shouldSelectChains() {
        Chain chain1 = Selection.on(protein)
                .asChain();
        System.out.println(chain1);
    }
}
