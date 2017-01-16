package de.bioforscher.jstructure.reconstruction.sidechain.pulchra;

import de.bioforscher.jstructure.mathematics.LinearAlgebraAtom;
import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.model.Combinatorics;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Element;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.family.AminoAcidFamily;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import de.bioforscher.jstructure.reconstruction.ReconstructionAlgorithm;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.*;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * Reduced implementation of the Pulchra algorithm.<br />
 * Intellectual and implementation credit to:<br />
 * <b>[Rotkiewicz, 2008] - Fast procedure for reconstruction of full-atom protein models from reduced representations, 10.1002/jcc.20906</b><br /><br />
 *
 * Places side chains by finding suitable rotamers from a library.<br /><br />
 *
 * A rather rudimentary adaptation of the PULCHRA algorithm for side getChain
 * placement.<br />
 * Actually, several other features are available in the native impl such as
 * backbone placement or several refinement steps. Long term, they may be added
 * to this implementation. For now only plain side getChain reconstruction is
 * supported.<br />
 * <br />
 *
 * original comment:
 *
 * <pre>
 * PULCHRA Protein Chain Restoration Algorithm
 *
 * Version 3.04 July 2007 Contact: Piotr Rotkiewicz, piotr -at- pirx -dot- com
 * </pre>
 *
 * J Comput Chem. 2008 July 15; 29(9): 1460–1465. doi:10.1002/jcc.20906.
 */
public class PULCHRA implements ReconstructionAlgorithm {
    private static final double BIN_SIZE = 0.3;
    private static final String BASE_PATH = "pulchra/";
    /*
     * side-getChain library
     */
    private static final String ROT_STAT_IDX_LIBRARY = BASE_PATH + "rot_data_idx.dat";
    private static final String ROT_STAT_COORDS_LIBRARY = BASE_PATH + "rot_data_coords.dat";
    /*
     * backbone library
     */
    private static final String NCO_LIBRARY = BASE_PATH + "nco_data.dat";
    private static final String NCO_PRO_LIBRARY = BASE_PATH + "nco_data_pro.dat";

    private static List<int[]> rotStatIdx;
    private static List<double[]> rotStatCoords;
    private static Map<int[], double[]> ncoStat;
    private static Map<int[], double[]> ncoStatPro;
    private Map<int[], List<double[]>> ncoLibrary;
    private Map<int[], List<double[]>> ncoProLibrary;

    public PULCHRA() {
        if(rotStatIdx == null) {
            initializeLibrary();
        }
    }

    @Override
    public void reconstruct(Protein protein) {
        reconstructBackbone(protein);
        reconstructSidechains(protein);
    }

    private void reconstructBackbone(Protein protein) {

    }

    private void reconstructSidechains(Protein protein) {
        // TODO again reconstruct missing leading/tailing residues
        Combinatorics.fragmentsOf(Selection.on(protein)
                .aminoAcids()
                .asFilteredGroups()
                .collect(Collectors.toList()), 4)
                .forEach(fragment -> {
                    Group residueToReconstruct = fragment.getElement(2);
                    AminoAcidFamily aminoAcid = residueToReconstruct.getGroupInformation().getAminoAcidFamily();
                    if (aminoAcid.equals(AminoAcidFamily.GLYCINE) || aminoAcid.equals(AminoAcidFamily.UNKNOWN)) {
                        return;
                    }

                    double[] ca_p2 = Selection.on(fragment.getElement(0))
                            .alphaCarbonAtoms()
                            .asAtom()
                            .getCoordinates();
                    double[] ca_p1 = Selection.on(fragment.getElement(1))
                            .alphaCarbonAtoms()
                            .asAtom()
                            .getCoordinates();
                    double[] ca_tr = Selection.on(residueToReconstruct)
                            .alphaCarbonAtoms()
                            .asAtom()
                            .getCoordinates();
                    double[] ca_n1 = Selection.on(fragment.getElement(3))
                            .alphaCarbonAtoms()
                            .asAtom()
                            .getCoordinates();

                    double d13 = LinearAlgebra3D.distance(ca_p2, ca_tr);
                    double d24 = LinearAlgebra3D.distance(ca_p1, ca_n1);
                    double d14 = LinearAlgebra3D.distance14(ca_p2, ca_p1, ca_tr, ca_n1);
                    int[] residueBinning = binResidues(d13, d24, d14);
//                    System.out.println("residueBinning : " + Arrays.toString(residueBinning));

                    // find closest rotamer conformation
                    int[] bestMatchingRotamer = null;
                    double bestMatchingRotamerDistance = Double.MAX_VALUE;
                    int aminoAcidIndex = getAminoAcidIndex(residueToReconstruct);
//                    System.out.println("aminoAcidIndex : " + aminoAcidIndex);

                    for (int[] rotamer : rotStatIdx) {
//                        System.out.println("rotamer : " + Arrays.toString(rotamer));
                        if (rotamer[0] != aminoAcidIndex) {
                            continue;
                        }

                        double rotamerDistance = difference(rotamer, residueBinning);
//                        System.out.println("distance : " + rotamerDistance);
                        if (rotamerDistance < bestMatchingRotamerDistance) {
                            bestMatchingRotamerDistance = rotamerDistance;
                            bestMatchingRotamer = rotamer;
                        }
                    }

                    // new rebuild
                    double[][] rotation = rotation(ca_p1, ca_tr, ca_n1);

                    int pos = Objects.requireNonNull(bestMatchingRotamer)[5];
                    int nsc = (int) aminoAcid.sideChainAtomNames().count();

                    // all atoms within the coordinate file describing this getResidue
                    for (int i = 0; i < nsc; i++) {
                        String atomName = aminoAcid.getSideChainAtomNames().get(i);
                        //TODO as pdbSerial cannot be decided yet, we need to update them later or externally, boilerplately
                        Atom reconstructedAtom = new Atom(atomName, 0,
                                Element.valueOfIgnoreCase(atomName.substring(0, 1)), rotStatCoords.get(pos + i + 1));
                        // transform atom
                        reconstructedAtom = new LinearAlgebraAtom.Transformation(ca_tr,
                                rotation).transformCoordinates(reconstructedAtom);
                        residueToReconstruct.addAtom(reconstructedAtom);
                    }
                });
    }

    private double[][] rotation(double[] ca_p1, double[] ca_tr, double[] ca_n1) {
        double[] difference21 = LinearAlgebra3D.normalize(LinearAlgebra3D.subtract(ca_tr, ca_p1));
        double[] difference23 = LinearAlgebra3D.normalize(LinearAlgebra3D.subtract(ca_tr, ca_n1));
        double[] difference13m = LinearAlgebra3D.normalize(LinearAlgebra3D.subtract(difference21, difference23));
        double[] difference13p = LinearAlgebra3D.normalize(LinearAlgebra3D.add(difference21, difference23));

        return new double[][] { difference13m, difference13p,
                { difference13m[1] * difference13p[2] - difference13m[2] * difference13p[1],
                        difference13m[2] * difference13p[0] - difference13m[0] * difference13p[2],
                        difference13m[0] * difference13p[1] - difference13m[1] * difference13p[0]
                }};
    }

    private double difference(int[] v1, int[] v2) {
        return Math.abs(v1[0] - v2[0]) + Math.abs(v1[1] - v2[1]) + 0.2 * Math.abs(v1[2] - v2[2]);
    }

    private int getAminoAcidIndex(Group residue) {
        final String aminoAcids = "GASCVTIPMDNLKEQRHFYWX";
        return aminoAcids.indexOf(residue.getGroupInformation().getOneLetterCode());
    }

    private int[] binResidues(double d13_1, double d13_2, double d14) {
        int bin13_1 = (int) ((d13_1 - 4.6) / BIN_SIZE);
        int bin13_2 = (int) ((d13_2 - 4.6) / BIN_SIZE);
        int bin14 = (int) ((d14 + 11.0) / BIN_SIZE);

        bin13_1 = LinearAlgebra3D.capToInterval(0, bin13_1, 9);
        bin13_2 = LinearAlgebra3D.capToInterval(0, bin13_2, 9);
        bin14 = LinearAlgebra3D.capToInterval(0, bin14, 73);

        return new int[] { bin13_1, bin13_2, bin14 };
    }

    private static final Predicate<String> LINE_FILTER = line -> line.contains("},");
    private static final Function<String, String[]> STRING_MAPPER = line -> line.replace("{", "").replace("}",
            "").split(",");

    private synchronized void initializeLibrary() {
        //TODO this could use some sophisticated parser - further versions of the data could easily break the code

        // parse sidechain indices
        InputStream idxIs = getResourceAsStream(ROT_STAT_IDX_LIBRARY);
        rotStatIdx = new BufferedReader(new InputStreamReader(idxIs)).lines()
            // select lines which do not contain information
            .filter(LINE_FILTER)
            // remove padding
            .map(STRING_MAPPER)
            .map(this::parseRotStatIdxLine)
            .collect(Collectors.toList());

        // parse side getChain coordinates
        InputStream coordsIs = getResourceAsStream(ROT_STAT_COORDS_LIBRARY);
        rotStatCoords = new BufferedReader(new InputStreamReader(coordsIs)).lines()
            // select lines which do not contain information
            .filter(LINE_FILTER)
            // remove padding
            .map(STRING_MAPPER)
            .map(this::parseRotStatCoordsLine)
            .collect(Collectors.toList());

        // parse backbone library for proline
        InputStream ncoData = getResourceAsStream(NCO_LIBRARY);
        ncoLibrary = new BufferedReader(new InputStreamReader(ncoData)).lines()
                .collect(NCOConsumer::new, NCOConsumer::accept, NCOConsumer::combine)
                .getData();

        // parse backbone lirbary for non-prolines
        InputStream ncoProDat = getResourceAsStream(NCO_PRO_LIBRARY);
        ncoProLibrary = new BufferedReader(new InputStreamReader(ncoProDat)).lines()
                .collect(NCOConsumer::new, NCOConsumer::accept, NCOConsumer::combine)
                .getData();
    }

    static class NCOConsumer implements Consumer<String> {
        private Map<int[], List<double[]>> data;
        private int[] currentBinIndex;
        private List<double[]> currentPieceOfData;

        NCOConsumer() {
            this.data = new HashMap<>();
        }

        @Override
        public void accept(String line) {
            // entry starting
            if(line.contains("{")) {
                currentBinIndex = splitToBinIndex(line);
                currentPieceOfData = new ArrayList<>();
            } else if(line.contains("}}")) {
                // entry ending
                data.put(currentBinIndex, currentPieceOfData);
            } else {
                // coordinate line
                currentPieceOfData.add(splitToData(line));
            }
        }

        private double[] splitToData(String line) {
            return Pattern.compile(",").splitAsStream(line)
                                       .map(String::trim)
                                       .mapToDouble(Double::valueOf)
                                       .toArray();
        }

        private int[] splitToBinIndex(String line) {
            return Pattern.compile(",").splitAsStream(line.split("\\{  \\{")[1].split("},")[0])
                                       .map(String::trim)
                                       .mapToInt(Integer::valueOf)
                                       .toArray();
        }

        void combine(NCOConsumer other) {
            data.putAll(other.data);
        }

        /**
         * Gets the composed data as map.
         * @return the composed data
         */
        Map<int[], List<double[]>> getData() {
            return data;
        }
    }

    private double[] parseRotStatCoordsLine(String[] tmp) {
        return Arrays.stream(tmp)
                     .map(String::trim)
                     .filter(line -> !line.isEmpty())
                     .mapToDouble(Double::valueOf)
                     .toArray();
    }

    private int[] parseRotStatIdxLine(String[] tmp) {
        return Arrays.stream(tmp)
                     .map(String::trim)
                     .filter(line -> !line.isEmpty())
                     .mapToInt(Integer::valueOf)
                     .toArray();
    }

    private InputStream getResourceAsStream(String filepath) {
        return Objects.requireNonNull(Thread.currentThread().getContextClassLoader().getResourceAsStream(filepath),
                "failed to findAny resource as InputStream");
    }
}
