package studies.gmlvq.hds;

import studies.StudyConstants;
import studies.gmlvq.AbstractItemSetMinerDataSetComposer;

import java.io.IOException;

/**
 * Create the new data set for the HDS-case.
 * Created by bittrich on 5/22/17.
 */
class S01_ComposeItemSetMinerDataSet extends AbstractItemSetMinerDataSetComposer {
    private static final String FILEPATH = StudyConstants.GMLVQ_MAIN + "data/itemset_miner/PF00089/";
    private static final String POSITIVE_FILENAME = FILEPATH + "positives.txt";
    private static final String NEGATIVE_FILENAME = FILEPATH + "negatives.txt";
    private static final String OUTPUT = FILEPATH + "PF00089.arff";

    public static void main(String[] args) throws IOException {
        new S01_ComposeItemSetMinerDataSet(POSITIVE_FILENAME, NEGATIVE_FILENAME, OUTPUT);
    }

    private S01_ComposeItemSetMinerDataSet(String positiveFilename, String negativeFilename, String outputFilename) throws IOException {
        super(positiveFilename, negativeFilename, outputFilename);
    }
}
