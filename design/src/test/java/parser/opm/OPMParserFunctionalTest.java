package parser.opm;

import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import design.DesignConstants;
import design.ProteinSource;
import design.parser.opm.OPMParser;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;
import java.nio.file.Paths;

/**
 * Ensures integrity of the OPM parser.
 * Created by S on 29.10.2016.
 */
public class OPMParserFunctionalTest {
    private Protein protein1xio;
    private String pdbId;

    @Before
    public void setup() {
        pdbId = "1xio";
        protein1xio = ProteinParser.parseProteinById(pdbId);
    }

    @Test
    public void shouldParseOPMFile() throws IOException {
        OPMParser.parse(protein1xio, Paths.get(DesignConstants.OPM_RAW_DIR + pdbId + DesignConstants.OPM_SUFFIX));
    }

    @Test
    public void shouldParseAllOPMFiles() throws IOException {
        ProteinSource.loadProteins(false, false, true);
    }
}
