package de.bioforscher.jstructure.mathematics.mds;

import de.bioforscher.testutil.TestUtils;
import org.junit.Assert;
import org.junit.Test;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Tests the MDS implementation.
 * Created by S on 27.09.2016.
 */
public class MultiDimensionalScalingTest {
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
            double[][] dataPoints = parseDistancesFromPocketAlignCSV(TestUtils.getResourceAsLines(resource));
            List<String> dataLabels = parseLabels(TestUtils.getResourceAsLines(resource));
            List<double[]> embeddedDataPoints = mds.computeEmbedding(dataPoints, 3);
            for(int i = 0; i < dataLabels.size(); i++) {
                String s = Arrays.toString(embeddedDataPoints.get(i));
                System.out.println(dataLabels.get(i) + "," + s.replace("[", "").replace("]", "").replace(" ", ""));
            }
        } catch (NullPointerException e) {
            Assert.fail("could not find test resource: " + e.getLocalizedMessage());
        }
    }

    private List<String> parseLabels(List<String> lines) {
        String line = lines.get(0);
        return Arrays.stream(line.split(",")).collect(Collectors.toList());
    }

    private double[][] parseDistancesFromPocketAlignCSV(List<String> lines) {
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