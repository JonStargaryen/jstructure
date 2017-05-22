package studies.gmlvq.cmh;

import studies.StudyConstants;
import studies.gmlvq.AbstractItemSetMinerDataSetComposer;

import java.io.IOException;

/**
 * Create some data set for the GMLVQ_MAIN-project. Handles the metal-binding motif.
 * Created by bittrich on 4/11/17.
 */
class S01_ComposeItemSetMinerDataSet extends AbstractItemSetMinerDataSetComposer {
    private static final String FILEPATH = StudyConstants.GMLVQ_MAIN + "data/itemset_miner/PF00127/";
    private static final String POSITIVE_FILENAME = FILEPATH + "positives.txt";
    private static final String NEGATIVE_FILENAME = FILEPATH + "negatives.txt";

    public static void main(String[] args) throws IOException {
        new S01_ComposeItemSetMinerDataSet(POSITIVE_FILENAME, NEGATIVE_FILENAME);
    }

    private S01_ComposeItemSetMinerDataSet(String positiveFilename, String negativeFilename) throws IOException {
        super(positiveFilename, negativeFilename);
    }
}
