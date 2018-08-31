package de.bioforscher.jstructure.efr.analysis.aars;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.StringJoiner;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * Compute the EF conservation of each column in the MSA. Are some positions highly EFing or are the predictions
 * randomly distributed along the sequence?
 * TODO: a possible extension would be to fragmentize the dataset further by aaRS type
 * best decision thresholds: 0.163 for the estimated probabilities.
 */
public class A02_EFConservationPerColumn {
    private static final double CUTOFF = 0.163;
    private static List<Column> class1Columns;
    private static List<Column> class2Columns;

    public static void main(String[] args) throws IOException {
        // load column conservation
        class1Columns = extractSequenceConservationFromFile("/home/bittrich/git/aars_data/T05_msa/C1_conservation.txt").stream()
                .map(Column::new)
                .collect(Collectors.toList());
        class2Columns = extractSequenceConservationFromFile("/home/bittrich/git/aars_data/T05_msa/C2_conservation.txt").stream()
                .map(Column::new)
                .collect(Collectors.toList());

        Files.list(Paths.get("/home/bittrich/git/phd_sb_repo/data/aars/EFoldMine_code/ef-predictions/"))
                .forEach(predictionPath -> {
                    try {
                        String aarsId = predictionPath.toFile().getName().split("\\.")[0];
                        System.out.println("processing " + aarsId);
//                        String[] split = aarsId.split("_");
//                        String pdbId = split[0];
//                        String chainId = split[1];
                        String className = Files.exists(Paths.get("/home/bittrich/git/aars_data/T06_renumbered_structures/C1/" + aarsId + ".pdb")) ? "C1" : "C2";
//                        Path originalStructurePath = Paths.get("/home/bittrich/git/aars_data/T06_renumbered_structures/" + className + "/" + aarsId + ".pdb");
                        Path renumberedStructurePath = Paths.get("/home/bittrich/git/aars_data/T06_renumbered_structures/" + className + "/renum/" + aarsId + "_renum.pdb");
//                        Structure originalStructure = StructureParser.fromPath(originalStructurePath).parse();
                        Structure renumberedStructure = StructureParser.fromPath(renumberedStructurePath).parse();

                        List<Column> columns = className.equals("C1") ? class1Columns : class2Columns;

                        List<Double> values = Files.lines(predictionPath)
                                .skip(1)
                                .filter(line -> !line.startsWith("*"))
                                .map(line -> line.split("\\s+")[1])
                                .map(Double::valueOf)
                                .collect(Collectors.toList());
//                        String sequence = Files.lines(predictionPath)
//                                .skip(1)
//                                .filter(line -> !line.startsWith("*"))
//                                .map(line -> line.split("\\s+")[0])
//                                .collect(Collectors.joining());

                        // assign EFpredictions to correct column
                        for(int i = 0; i < values.size(); i++) {
                            double value = values.get(i);
//                            AminoAcid originalAminoAcid = originalStructure.getFirstChain().getAminoAcids().get(i);
                            AminoAcid renumberedAminoAcid = renumberedStructure.getFirstChain().getAminoAcids().get(i);

                            // check for sanity
//                            if(!originalAminoAcid.getThreeLetterCode().equals(renumberedAminoAcid.getThreeLetterCode())) {
//                                System.out.println(originalAminoAcid.toString() + " -> " + renumberedAminoAcid.toString() + " : " + value);
//                            }

                            Column column = columns.get(renumberedAminoAcid.getResidueIdentifier().getResidueNumber() - 1);
                            column.addEfPrediction(value);
                        }

                        // assign functional annotation
                        // from individual binding sites (this will miss structures where no ligands are present)
//                        Path contactPath = Paths.get("/home/bittrich/git/aars_data/T09_interactions/" + className + "/contacts_per_structure.txt");
//                        List<Integer> functionalPositions = Files.lines(contactPath)
//                                .filter(line -> line.startsWith(aarsId))
//                                .map(line -> line.split("\t")[1])
//                                .map(line -> line.substring(1, line.length() - 1))
//                                .flatMap(line -> Pattern.compile(",").splitAsStream(line))
//                                .map(Integer::valueOf)
//                                .collect(Collectors.toList());
                        // consider all every observed binding sites across all aaRS types
                        Path contactPath = Paths.get("/home/bittrich/git/aars_data/T09_interactions/" + className + "/contacts_all.txt");
                        List<Integer> functionalPositions = Files.lines(contactPath)
                                .map(line -> line.split("\t")[1])
                                .map(line -> line.substring(1, line.length() - 1))
                                .flatMap(line -> Pattern.compile(",").splitAsStream(line))
                                .map(Integer::valueOf)
                                .collect(Collectors.toList());

                        for(int functionalPosition : functionalPositions) {
                            columns.get(functionalPosition - 1).setFunctional(true);
                        }
                    } catch (IOException e) {
                        throw new UncheckedIOException(e);
                    }
                });

        System.out.println("class I");
        String class1Line = composeClassSpecificOutput(class1Columns,
                Paths.get("/home/bittrich/git/phd_sb_repo/data/aars/EFoldMine_code/C1_ef_frequencies.csv"),
                81);
        System.out.println("class II");
        String class2Line = composeClassSpecificOutput(class2Columns,
                Paths.get("/home/bittrich/git/phd_sb_repo/data/aars/EFoldMine_code/C2_ef_frequencies.csv"),
                75);

        Files.write(Paths.get("/home/bittrich/git/phd_sb_repo/data/aars/EFoldMine_code/summary.csv"),
                ("class;total;func;ef_cons;ef_avg;intersect_cons;intersect_avg" + System.lineSeparator() +
                "class1;" + class1Line + System.lineSeparator() +
                "class2;" + class2Line).getBytes()
        );
    }

    private static String composeClassSpecificOutput(List<Column> columns,
                                                     Path outputPath,
                                                     int numberOfRepresentatives) {
        StringJoiner output = new StringJoiner(System.lineSeparator());
        output.add("resNum;func;ef_cons;ef_avg;conservation;n_data;ef_score_avg;ef_freq");
        for(int i = 0; i < columns.size(); i++) {
            Column column = columns.get(i);
            String line = (i + 1) + ";" +
                    column.isFunctional() + ";" +
                    column.isConservedEf() + ";" +
                    column.isOnAverageEf() + ";" +
                    StandardFormat.format(column.getConservationValue()) + ";" +
                    column.getNumberOfEfData() + ";" +
                    StandardFormat.format(column.getAverageEfPrediction()) + ";" +
                    StandardFormat.format(column.getNumberOfEfPredictions() / (double) numberOfRepresentatives);
            System.out.println(line);
            output.add(line);
        }

        // determine class-specific summary - only consider columns in the protozyme
        boolean class1 = columns.size() != 2515;
        List<Column> filteredColumns = new ArrayList<>();
        for(int i = 0; i < columns.size(); i++) {
            Column column = columns.get(i);
            if((class1 && (i > 254 && i < 337)) || (!class1 && (i > 647 && i < 718))) {
                filteredColumns.add(column);
            }
        }
        int total = filteredColumns.size();
        int functional = (int) filteredColumns.stream()
                .filter(Column::isFunctional)
                .count();
        int ef_cons = (int) filteredColumns.stream()
                .filter(Column::isConservedEf)
                .count();
        int ef_avg = (int) filteredColumns.stream()
                .filter(Column::isOnAverageEf)
                .count();
        int intersect_cons = (int) filteredColumns.stream()
                .filter(column -> column.isFunctional() && column.isConservedEf())
                .count();
        int intersect_avg = (int) filteredColumns.stream()
                .filter(column -> column.isFunctional() && column.isOnAverageEf())
                .count();

        try {
            Files.write(outputPath, output.toString().getBytes());
            return total + ";" + functional + ";" + ef_cons + ";" + ef_avg + ";" + intersect_cons + ";" + intersect_avg;
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    private static List<Double> extractSequenceConservationFromFile(String filename) throws IOException {
        String line = Files.lines(Paths.get(filename))
                .filter(l -> l.startsWith("BAR_GRAPH\tConservation\tConservation of total alignment less than 25% gaps\t"))
                .findFirst()
                .get()
                .split("\t")[3];

        return Pattern.compile("\\|").splitAsStream(line)
                .map(entry -> entry.split(",")[0])
                .map(Double::valueOf)
                .collect(Collectors.toList());
    }

    static class Column {
        private double conservationValue;
        private List<Double> efPredictions;
        private boolean functional;

        Column(double conservationValue) {
            this.conservationValue = conservationValue;
            this.efPredictions = new ArrayList<>();
        }

        public double getConservationValue() {
            return conservationValue;
        }

        public List<Double> getEfPredictions() {
            return efPredictions;
        }

        public int getNumberOfEfData() {
            return efPredictions.size();
        }

        public double getAverageEfPrediction() {
            return efPredictions.stream()
                    .mapToDouble(Double::valueOf)
                    .average()
                    .orElse(0.0);
        }

        public boolean isOnAverageEf() {
            return getAverageEfPrediction() > CUTOFF;
        }

        public boolean isConservedEf() {
            return isOnAverageEf() && getNumberOfEfData() > 4;
        }

        public double getNumberOfEfPredictions() {
            return efPredictions.stream()
                    .filter(efPrediction -> efPrediction > CUTOFF)
                    .count();
        }

        public void addEfPrediction(double value) {
            this.efPredictions.add(value);
        }

        public boolean isFunctional() {
            return functional;
        }

        public void setFunctional(boolean functional) {
            this.functional = functional;
        }
    }
}
