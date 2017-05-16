package bioforscher.kinks;

import bioforscher.Constants;

import java.nio.file.Paths;
import java.util.StringJoiner;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Compute several features and compose arff.
 * Created by bittrich on 2/14/17.
 */
@Deprecated
public class S02_ComposeKinkFinderArff {
    public static void main(String[] args) {
        StringJoiner header = new StringJoiner(System.lineSeparator());
        /*
         * @RELATION iris
         * @ATTRIBUTE sepallength  NUMERIC
         * @ATTRIBUTE sepalwidth   NUMERIC
         * @ATTRIBUTE petallength  NUMERIC
         * @ATTRIBUTE petalwidth   NUMERIC
         * @ATTRIBUTE class        {Iris-setosa,Iris-versicolor,Iris-virginica}
         */
        header.add("@RELATION kinks");
        header.add("@ATTRIBUTE pdbId     STRING");
        header.add("@ATTRIBUTE chainId   STRING");
        header.add("@ATTRIBUTE start     NUMERIC");
        header.add("@ATTRIBUTE end       NUMERIC");
        header.add("@ATTRIBUTE kink      NUMERIC");
        header.add("@ATTRIBUTE gly       NUMERIC");
        header.add("@ATTRIBUTE pro       NUMERIC");
        header.add("@ATTRIBUTE rasaL4    NUMERIC");
        header.add("@ATTRIBUTE rasaL3    NUMERIC");
        header.add("@ATTRIBUTE rasaL2    NUMERIC");
        header.add("@ATTRIBUTE rasaL1    NUMERIC");
        header.add("@ATTRIBUTE rasa      NUMERIC");
        header.add("@ATTRIBUTE rasaR1    NUMERIC");
        header.add("@ATTRIBUTE rasaR2    NUMERIC");
        header.add("@ATTRIBUTE rasaR3    NUMERIC");
        header.add("@ATTRIBUTE rasaR4    NUMERIC");
        header.add("@ATTRIBUTE plipL4    NUMERIC");
        header.add("@ATTRIBUTE plipL3    NUMERIC");
        header.add("@ATTRIBUTE plipL2    NUMERIC");
        header.add("@ATTRIBUTE plipL1    NUMERIC");
        header.add("@ATTRIBUTE plip      NUMERIC");
        header.add("@ATTRIBUTE plipR1    NUMERIC");
        header.add("@ATTRIBUTE plipR2    NUMERIC");
        header.add("@ATTRIBUTE plipR3    NUMERIC");
        header.add("@ATTRIBUTE plipR4    NUMERIC");
        header.add("@ATTRIBUTE hbondL4    NUMERIC");
        header.add("@ATTRIBUTE hbondL3    NUMERIC");
        header.add("@ATTRIBUTE hbondL2    NUMERIC");
        header.add("@ATTRIBUTE hbondL1    NUMERIC");
        header.add("@ATTRIBUTE hbond      NUMERIC");
        header.add("@ATTRIBUTE hbondR1    NUMERIC");
        header.add("@ATTRIBUTE hbondR2    NUMERIC");
        header.add("@ATTRIBUTE hbondR3    NUMERIC");
        header.add("@ATTRIBUTE hbondR4    NUMERIC");
        header.add("@ATTRIBUTE sideL4    NUMERIC");
        header.add("@ATTRIBUTE sideL3    NUMERIC");
        header.add("@ATTRIBUTE sideL2    NUMERIC");
        header.add("@ATTRIBUTE sideL1    NUMERIC");
        header.add("@ATTRIBUTE side      NUMERIC");
        header.add("@ATTRIBUTE sideR1    NUMERIC");
        header.add("@ATTRIBUTE sideR2    NUMERIC");
        header.add("@ATTRIBUTE sideR3    NUMERIC");
        header.add("@ATTRIBUTE sideR4    NUMERIC");
        header.add("@ATTRIBUTE motifs    STRING");
        header.add("@ATTRIBUTE error     NUMERIC");
        header.add("@ATTRIBUTE angle     NUMERIC");

        /*
         * @DATA
         * 5.1,3.5,1.4,0.2,Iris-setosa
         * 4.9,3.0,1.4,0.2,Iris-setosa
         * 4.7,3.2,1.3,0.2,Iris-setosa
         * 4.6,3.1,1.5,0.2,Iris-setosa
         * 5.0,3.6,1.4,0.2,Iris-setosa
         * 5.4,3.9,1.7,0.4,Iris-setosa
         * 4.6,3.4,1.4,0.3,Iris-setosa
         * 5.0,3.4,1.5,0.2,Iris-setosa
         * 4.4,2.9,1.4,0.2,Iris-setosa
         * 4.9,3.1,1.5,0.1,Iris-setosa
         */
        header.add("@DATA");
        String output = Constants.lines(Paths.get(Constants.STATISTICS_PATH + "kink_finder_results_all.csv"))
                .filter(line -> !line.startsWith("pdbId"))
                .map(line -> line.replace(",", "|"))
                .map(line -> line.split("\t"))
                .map(split -> Stream.of(split).collect(Collectors.joining(",")))
                .peek(System.out::println)
                .collect(Collectors.joining(System.lineSeparator(), header.toString() + System.lineSeparator(), ""));

        Constants.write(Paths.get(Constants.STATISTICS_PATH + "kink_finder_results_all.arff"), output.getBytes());
    }
}
