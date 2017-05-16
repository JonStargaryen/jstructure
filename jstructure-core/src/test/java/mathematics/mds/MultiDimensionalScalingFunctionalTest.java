package mathematics.mds;

import de.bioforscher.jstructure.mathematics.mds.MultiDimensionalScaling;
import org.junit.Assert;
import org.junit.Test;
import util.TestUtils;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by S on 27.09.2016.
 */
public class MultiDimensionalScalingFunctionalTest {
    private static final String POCKET_ALIGN_CSV_1_PATH = "mathematics/mds/scheme1.csv";
    private static final String POCKET_ALIGN_CSV_2_PATH = "mathematics/mds/scheme2.csv";

    @Test
    public void shouldRunMultiDimensionalScaling() {
        runMultiDimensionalScalingForResource(POCKET_ALIGN_CSV_1_PATH);
        runMultiDimensionalScalingForResource(POCKET_ALIGN_CSV_2_PATH);
    }

    private void runMultiDimensionalScalingForResource(String resource) {
        System.out.println("processing " + resource);
        try {
            MultiDimensionalScaling mds = new MultiDimensionalScaling();
            double[][] dataPoints = parseDistancesFromPocketAlignCSV(TestUtils.getResourceAsFilepath(resource));
            List<String> dataLabels = parseLabels(TestUtils.getResourceAsFilepath(resource));
            List<double[]> embeddedDataPoints = mds.computeEmbedding(dataPoints, 3);
            for(int i = 0; i < dataLabels.size(); i++) {
                String s = Arrays.toString(embeddedDataPoints.get(i));
                System.out.println(dataLabels.get(i) + "," + s.replace("[", "").replace("]", "").replace(" ", ""));
            }
        } catch (IOException e) {
            Assert.fail("failed with IOException: " + e.getLocalizedMessage());
        } catch (NullPointerException e) {
            Assert.fail("could not find test resource: " + e.getLocalizedMessage());
        }
    }

    private List<String> parseLabels(String filepath) throws IOException {
        String line = Files.readAllLines(new File(filepath).toPath()).get(0);
        return Arrays.stream(line.split(",")).collect(Collectors.toList());
    }

    private double[][] parseDistancesFromPocketAlignCSV(String filepath) throws IOException {
        List<String> lines = Files.readAllLines(new File(filepath).toPath());
        double[][] data = new double[lines.size() - 1][];

        // skip first line
        for(int i = 0; i < lines.size(); i++) {
            if(i == 0) {
                continue;
            }
            String[] tmpLine = lines.get(i).split(",");
            double[] tmpData = new double[lines.size() - 1];
            for(int j = 0; j < lines.size(); j++) {
                if(j == 0) {
                    continue;
                }
                if(i == j) {
                    tmpData[j - 1] = 0.0;
                    continue;
                }
                tmpData[j - 1] = Double.valueOf(tmpLine[j]);
            }
            data[i - 1] = tmpData;
        }

        return data;
    }
}