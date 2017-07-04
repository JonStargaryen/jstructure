package studies.gmlvq.fingerprint;

import studies.StudyConstants;
import weka.classifiers.functions.GMLVQ;
import weka.core.Instances;
import weka.filters.Filter;
import weka.filters.unsupervised.attribute.Remove;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.stream.IntStream;

/**
 * Train a model on a group of data sets using the Java API.
 * Created by bittrich on 7/4/17.
 */
public class S02_TrainGmlvqForFingerprintData {
    private static final int NUMBER_OF_PLIP_FEATURES = 10;

    public static void main(String[] args) throws IOException {
        // traverse directory with arff files
        Files.list(StudyConstants.FINGERPRINT_MINER_PFAM_ARFF)
                .map(Path::toFile)
                .forEach(S02_TrainGmlvqForFingerprintData::handlePath);
    }

    private static void handlePath(File file) {
        System.out.println(file);
        try {
            BufferedReader reader = new BufferedReader(new FileReader(file));
            // load data into weka's data structure
            Instances data = new Instances(reader);
            reader.close();

            // data sets can have varying number of features - determine the number how often they repeat
            int numberOfResidues = (data.numAttributes() - 2) / NUMBER_OF_PLIP_FEATURES;
            // determine indices of attributes to remove
            int[] indicesToRemove = IntStream.concat(IntStream.of(0),
                    IntStream.range(0, numberOfResidues)
                            .map(residueNumber -> residueNumber * NUMBER_OF_PLIP_FEATURES)
                            .flatMap(residueNumber -> IntStream.of(4, 6, 7, 8, 9, 10).map(i -> residueNumber + i)))
                    .toArray();
            // remove ids
            Remove remove = new Remove();
            remove.setAttributeIndicesArray(indicesToRemove);
            remove.setInputFormat(data);
            // employ filer
            data = Filter.useFilter(data, remove);

            // mark last attribute as class
            data.setClassIndex(data.numAttributes() - 1);

            // create and build GMLVQ instance
            GMLVQ gmlvq = new GMLVQ();
            gmlvq.setVisualization(false);
            gmlvq.buildClassifier(data);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
