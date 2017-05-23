package studies.gmlvq.cmh;

import studies.StudyConstants;
import studies.gmlvq.AbstractItemSetMinerDataSetComposer;

import java.io.IOException;

/**
 * Create some data set for the GMLVQ_MAIN-project. Handles the metal-binding motif.
 * Created by bittrich on 4/11/17.
 */
class S01_ComposeMetalDataSet extends AbstractItemSetMinerDataSetComposer {
    private static final String FILEPATH = StudyConstants.GMLVQ_MAIN + "data/itemset_miner/PF00127/";
    private static final String POSITIVE_FILENAME = FILEPATH + "positives.txt";
    private static final String NEGATIVE_FILENAME = FILEPATH + "negatives.txt";
    private static final String OUTPUT = FILEPATH + "PF00127.arff";

    public static void main(String[] args) throws IOException {
        new S01_ComposeMetalDataSet(POSITIVE_FILENAME, NEGATIVE_FILENAME, OUTPUT);
    }

    private S01_ComposeMetalDataSet(String positiveFilename, String negativeFilename, String outputFilename) throws IOException {
        super(positiveFilename, negativeFilename, outputFilename);
    }
}
