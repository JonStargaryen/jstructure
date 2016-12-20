package design;

import de.bioforscher.jstructure.model.structure.AminoAcid;

import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * A converter of minimal data to the <code>*.arff</code> format.
 * TODO Someday, a more generic implementation would be nice.
 * Created by S on 02.11.2016.
 */
public class ArffWriter {
    private static final String NOMINAL_AMINO_ACID_STRING = Arrays.stream(AminoAcid.values())
                                                                  .map(AminoAcid::getOneLetterCode)
                                                                  .collect(Collectors.joining(",", "{", "}"));

    /**
     * Exports the given data to a String. The last entry in the data argument will serve as class label.
     * @param name the name of the relation
     * @param data the data container
     * @return the String representation of the data
     */
    public static String convertToArff(String name, List<String> data) {
        final StringJoiner stringRepresentation = new StringJoiner(System.lineSeparator());

        /*
         * @RELATION iris
         * @ATTRIBUTE sepallength  NUMERIC
         * @ATTRIBUTE sepalwidth   NUMERIC
         * @ATTRIBUTE petallength  NUMERIC
         * @ATTRIBUTE petalwidth   NUMERIC
         * @ATTRIBUTE class        {Iris-setosa,Iris-versicolor,Iris-virginica}
         */
        stringRepresentation.add("@RELATION " + name);
        IntStream.range(1, data.get(0).split(",").length).forEach(i -> {
//            TODO formatting could be nicer
//            String attributeName = String.format("@ATTRIBUTE position%2d  " + NOMINAL_AMINO_ACID_STRING, i);
//            attributeName = attributeName.replace("position ", "position0");
//            stringRepresentation.add(attributeName);

            stringRepresentation.add(String.format("@ATTRIBUTE position%2d  NUMERIC", i).replace("position ", "position0"));
        });
//        stringRepresentation.add("@ATTRIBUTE class       {tm,ntm,trans}");
        stringRepresentation.add("@ATTRIBUTE class       {tm,ntm}");

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
        stringRepresentation.add("@DATA");
        data.stream()
            // TODO sometimes the string will be empty
            .filter(line -> !line.startsWith(","))
            .forEach(stringRepresentation::add);

        return stringRepresentation.toString();
    }

    /**
     * Exports the given data to the file system.
     * @param name the name of the relation
     * @param data the data container - each entry must provide a suitable {@link Object#toString()} implementation
     * @param path where to write?
     * @see #convertToArff(String, List)
     */
    public static void writeToArff(String name, List<String> data, Path path) {
        DesignConstants.write(path, convertToArff(name, data).getBytes());
    }
}
