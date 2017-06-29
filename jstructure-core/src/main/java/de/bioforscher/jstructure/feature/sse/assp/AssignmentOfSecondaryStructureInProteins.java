package de.bioforscher.jstructure.feature.sse.assp;

import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.feature.sse.SecondaryStructureElement;
import de.bioforscher.jstructure.mathematics.LinearAlgebra;
import de.bioforscher.jstructure.mathematics.TorsionAngles;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * An algorithm to determine the secondary structure of proteins by investigating the angles for 4 consecutive amino
 * acids.
 * @see <a href="http://nucleix.mbu.iisc.ac.in/assp/algorithm.html">http://nucleix.mbu.iisc.ac.in/assp/algorithm.html</a>
 * Created by bittrich on 6/28/17.
 */
@FeatureProvider(provides = { GenericSecondaryStructure.class, ASSPSecondaryStructure.class })
public class AssignmentOfSecondaryStructureInProteins extends AbstractFeatureProvider {
    private static final Logger logger = LoggerFactory.getLogger(AssignmentOfSecondaryStructureInProteins.class);

    @Override
    protected void processInternally(Protein protein) {
        protein.chainsWithAminoAcids()
                .forEach(this::processInternally);
    }

    private void processInternally(Chain chain) {
        chain.aminoAcids()
                .forEach(this::assignNeutralState);

        // create map of all amino acids of the given chain
        List<AminoAcid> aminoAcids = chain.aminoAcids()
                .collect(Collectors.toList());

        assignHelicalParameters(aminoAcids);
        assignBendAngles(aminoAcids);
        // assign helical characteristics (i.e. U U U) and obtain information on continuous stretches in the structure
        List<List<AminoAcid>> stretches = assignHelicalCharacteristics(aminoAcids);

        dropInitialAndTerminalCoil(stretches);
        assignFinalCharacteristic(stretches);
    }

    /**
     * Assigns the neutral state to an amino acid, i.e. marking its secondary structure as coil.
     * @param aminoAcid the residue to process
     */
    private void assignNeutralState(AminoAcid aminoAcid) {
        aminoAcid.getFeatureContainer().addFeature(new ASSPSecondaryStructure(this, SecondaryStructureElement.COIL));
    }

    private void assignHelicalParameters(List<AminoAcid> aminoAcids) {
        // for each fragment of 4 consecutive amino acids: evaluate their geometry to derive secondary structure
        for(int residueIndex = 0; residueIndex < aminoAcids.size() - 3; residueIndex++) {
            AminoAcid aa1 = aminoAcids.get(residueIndex);
            AminoAcid aa2 = aminoAcids.get(residueIndex + 1);
            AminoAcid aa3 = aminoAcids.get(residueIndex + 2);
            AminoAcid aa4 = aminoAcids.get(residueIndex + 3);

            if(aa1.getResidueNumber().getResidueNumber() + 1 != aa2.getResidueNumber().getResidueNumber() ||
                    aa2.getResidueNumber().getResidueNumber() + 1 != aa3.getResidueNumber().getResidueNumber() ||
                    aa3.getResidueNumber().getResidueNumber() + 1 != aa4.getResidueNumber().getResidueNumber()) {
                logger.warn("encountered non-consecutive amino acids {} {} {} {}",
                        aa1,
                        aa2,
                        aa3,
                        aa4);
                continue;
            }

            try {
                Atom ca1 = aa1.getCa();
                Atom ca2 = aa2.getCa();
                Atom ca3 = aa3.getCa();
                Atom ca4 = aa4.getCa();

                // bonds between atoms
                LinearAlgebra.PrimitiveDoubleArrayLinearAlgebra b1 = ca2.calculate().subtract(ca1);
                LinearAlgebra.PrimitiveDoubleArrayLinearAlgebra b2 = ca3.calculate().subtract(ca2);
                LinearAlgebra.PrimitiveDoubleArrayLinearAlgebra b3 = ca4.calculate().subtract(ca3);

                LinearAlgebra.PrimitiveDoubleArrayLinearAlgebra v1 = b1.subtract(b2);
                LinearAlgebra.PrimitiveDoubleArrayLinearAlgebra v2 = b2.subtract(b3);

                // direction cosine
                LinearAlgebra.PrimitiveDoubleArrayLinearAlgebra u = v1.vectorProduct(v2).normalize();

                // twist
                double rawTwist = v1.dotProduct(v2) / (v1.norm() * v2.norm());
                double costheta = 1 - rawTwist;
                // 57.2958 = 180 / Math.PI, used to convert to degree
                double t = Math.acos(rawTwist) * 57.2958;

                // virtual torsion angle vtor
                double vtor = TorsionAngles.torsionAngle(ca1, ca2, ca3, ca4);
                if(vtor < 0) {
                    vtor += 360;
                }

                // rise per residue
                double h = b2.dotProduct(u) / u.norm();
                if(h < 0) {
                    h = Math.abs(h);
                    t = 360 - t;
                }

                // radius
                double r = Math.sqrt(v1.norm() * v2.norm()) / (2 * costheta);

                ASSPSecondaryStructure secondaryStructure = getSecondaryStructure(aa1);
                secondaryStructure.setHelicalParameters(t, h, vtor, r, u.getValue());
            } catch (NullPointerException e) {
                logger.warn("missing alpha carbons on fragment {} {} {} {}",
                        aa1,
                        aa2,
                        aa3,
                        aa4);
            }
        }
    }

    private void assignBendAngles(List<AminoAcid> aminoAcids) {
        // for each fragment of 4 consecutive amino acids: evaluate their helical parameters to derive bend angle
        for(int residueIndex = 3; residueIndex < aminoAcids.size() - 3; residueIndex++) {
            AminoAcid aa4 = aminoAcids.get(residueIndex);
            ASSPSecondaryStructure secondaryStructure4 = getSecondaryStructure(aa4);

            AminoAcid aa1 = aminoAcids.get(residueIndex - 3);
            ASSPSecondaryStructure secondaryStructure1 = getSecondaryStructure(aa1);

            secondaryStructure4.setBa(LinearAlgebra.on(secondaryStructure1.getU()).angle(secondaryStructure4.getU()));
        }
    }

    private List<List<AminoAcid>> assignHelicalCharacteristics(List<AminoAcid> aminoAcids) {
        // found stretches
        List<List<AminoAcid>> stretches = new ArrayList<>();
        // currently in a continuous stretch
        List<AminoAcid> aminoAcidsInCurrentStretch = new ArrayList<>();
        int stretchId = 1;

        for(int residueIndex = 1; residueIndex < aminoAcids.size() - 1; residueIndex++) {
            if(isContinuousStretch(aminoAcids, residueIndex)) {
                // stretch criteria fulfilled
                AminoAcid aaPrevious = aminoAcids.get(residueIndex - 1);
                if(!aminoAcidsInCurrentStretch.contains(aaPrevious)) {
                    aminoAcidsInCurrentStretch.add(aaPrevious);
                }
                aminoAcidsInCurrentStretch.add(aminoAcids.get(residueIndex));
            } else {
                // stretch may just have ended
                if(!aminoAcidsInCurrentStretch.isEmpty()) {
                    // add terminal amino acid of stretch
                    System.out.println("stretch from " + aminoAcidsInCurrentStretch.get(1).getResidueNumber().getResidueNumber() + " to " +
                            (aminoAcidsInCurrentStretch.get(aminoAcidsInCurrentStretch.size() - 1).getResidueNumber().getResidueNumber() + 2));
                    // if so, determine secondary structure
                    for(int i = 0; i < aminoAcidsInCurrentStretch.size(); i++) {
                        ASSPSecondaryStructure secondaryStructure1 = getSecondaryStructure(aminoAcids.get(aminoAcids.indexOf(aminoAcidsInCurrentStretch.get(i))));
                        ASSPSecondaryStructure secondaryStructure2 = getSecondaryStructure(aminoAcids.get(aminoAcids.indexOf(aminoAcidsInCurrentStretch.get(i)) + 1));
                        double twist1 = secondaryStructure1.getT();
                        double twist2 = secondaryStructure2.getT();
                        double height1 = secondaryStructure1.getH();
                        double height2 = secondaryStructure2.getH();
                        double radius1 = secondaryStructure1.getR();
                        double radius2 = secondaryStructure2.getR();
                        double sumt = (twist1 + twist2) / 2;
                        double sumh = (height1 + height2) / 2;
                        double sumrad = (radius1 + radius2) / 2;
                        double twist_ss1 = 360 - twist1;

                        String[] helicalCharacteristics = { "U", "U", "U" };

                        if(i < aminoAcidsInCurrentStretch.size() - 2) {
                            ASSPSecondaryStructure secondaryStructure3 = getSecondaryStructure(aminoAcidsInCurrentStretch.get(i + 2));
                            double twist3 = secondaryStructure3.getT();
                            double sumt1 = (twist1 + twist2 + twist3) / 3;
                            double sumh1 = (height1 + height2 + secondaryStructure3.getH()) / 3;
                            double sumrad1 = (radius1 + radius2 + secondaryStructure3.getR()) / 3;
                            double sumt_l = 360 - (twist1 + twist2) / 2;
                            double sumt1_l = 360 - (twist1 + twist2 + twist3) / 3;

                            if(((((sumt > 93.6) && (sumt < 103.5)) || ((twist1 > 93.6) && (twist1 < 103.5)) || ((sumrad > 2.2) && (sumrad < 2.4)) || ((radius1 > 2.2) && (radius1 < 2.4))) && (((sumh > 0.9) && (sumh < 2.1)) || ((height1 > 0.9) && (height1 < 2.1)))) && (twist1 < 180)) {
                                helicalCharacteristics[0] = "A";
                            } else if(((((sumt_l > 93.6) && (sumt_l < 103.5)) || ((twist_ss1 > 93.6) && (twist_ss1 < 103.5)) || ((sumrad > 2.2) && (sumrad < 2.4)) || ((radius1 > 2.2) && (radius1 < 2.4))) && (((sumh > 0.9) && (sumh < 2.1)) || ((height1 > 0.9) && (height1 < 2.1)))) && (twist1 > 180)) {
                                helicalCharacteristics[0] = "a";
                            } else if((((radius1 > 1.1) && (radius1 < 1.9)) && ((twist1 > 223.6) && (twist1 < 253))) && ((height1 > 2.8) && (height1 <= 3.2))) {
                                helicalCharacteristics[0] = "P";
                            }

                            if(((((twist1 > 103.4) && (twist1 < 114.9)) || ((radius1 > 2) && (radius1 <= 2.2))) && ((height1 > 0.9) && (height1 < 2.1))) && (twist1 < 180)) {
                                helicalCharacteristics[1] = "G";
                            } else if(((((twist_ss1 > 103.4) && (twist_ss1 < 114.9)) || ((radius1 > 2) && (radius1 <= 2.2))) && ((height1 > 0.9) && (height1 < 2.1))) && (twist1 > 180)) {
                                helicalCharacteristics[1] = "g";
                            }

                            if(((((sumt1 > 77.9) && (sumt1 < 93.6)) || ((sumrad1 > 2.4) && (sumrad1 <= 2.7)) || ((twist1 > 77.9) && (twist1 < 93.6)) || ((radius1 > 2.4) && (radius1 <= 2.7))) && (((sumh1 > 0.9) && (sumh1 < 2.1)) || ((height1 > 0.9) && (height1 < 2.1)))) && (twist1 < 180)) {
                                helicalCharacteristics[2] = "I";
                            } else if(((((sumt1_l > 77.9) && (sumt1_l < 93.6)) || ((sumrad1 > 2.4) && (sumrad1 <= 2.7)) || ((twist_ss1 > 77.9) && (twist_ss1 < 93.6)) || ((radius1 > 2.4) && (radius1 <= 2.7))) && (((sumh1 > 0.9) && (sumh1 < 2.1)) || ((height1 > 0.9) && (height1 < 2.1)))) && (twist1 > 180)) {
                                helicalCharacteristics[2] = "i";
                            }
                        } else if(i == aminoAcidsInCurrentStretch.size() - 2) {
                            //TODO deviation for PHE-33: expected 'A U U', actually 'A U I'
                            double sumt_l = 360 - (twist1 + twist2) / 2;

                            if(((((sumt > 93.6) && (sumt < 103.5)) || ((twist1 > 93.6) && (twist1 < 103.5)) || ((sumrad > 2.2) && (sumrad < 2.4)) || ((radius1 > 2.2) && (radius1 < 2.4))) && (((sumh > 0.9) && (sumh < 2.1)) || ((height1 > 0.9) && (height1 < 2.1)))) && (twist1 < 180)) {
                                helicalCharacteristics[0] = "A";
                            } else if(((((sumt_l > 93.6) && (sumt_l < 103.5)) || ((twist_ss1 > 93.6) && (twist_ss1 < 103.5)) || ((sumrad > 2.2) && (sumrad < 2.4)) || ((radius1 > 2.2) && (radius1 < 2.4))) && (((sumh > 0.9) && (sumh < 2.1)) || ((height1 > 0.9) && (height1 < 2.1)))) && (twist1 > 180)) {
                                helicalCharacteristics[0] = "a";
                            } else if((((radius1 > 1.1) && (radius1 <= 1.9)) && ((twist1 > 223.6) && (twist1 < 253))) && ((height1 > 2.8) && (height1 <= 3.2))) {
                                helicalCharacteristics[0] = "P";
                            }

                            if(((((twist1 > 103.4) && (twist1 < 114.9)) || ((radius1 > 2) && (radius1 <= 2.2))) && ((height1 > 0.9) && (height1 < 2.1))) && (twist1 < 180)) {
                                helicalCharacteristics[1] = "G";
                            } else if(((((twist_ss1 > 103.4) && (twist_ss1 < 114.9)) || ((radius1 > 2) && (radius1 <= 2.2))) && ((height1 > 0.9) && (height1 < 2.1))) && (twist1 > 180)) {
                                helicalCharacteristics[1] = "g";
                            }

                            if(((((twist1 > 77.9) && (twist1 < 93.6)) || ((radius1 > 2.4) && (radius1 <= 2.7))) && (((height1 > 0.9) && (height1 < 2.1)))) && (twist1 < 180)) {
                                helicalCharacteristics[2] = "I";
                            } else if(((((twist_ss1 > 77.9) && (twist_ss1 < 93.6)) || ((radius1 > 2.4) && (radius1 <= 2.7))) && ((height1 > 0.9) && (height1 < 2.1))) && (twist1 > 180)) {
                                helicalCharacteristics[2] = "i";
                            }
                        } else {
                            if(((((twist1 > 93.6) && (twist1 < 103.4)) || ((radius1 > 2.2) && (radius1 < 2.4))) && ((height1 > 0.9) && (height1 < 2.1))) && (twist1 < 180)) {
                                helicalCharacteristics[0] = "A";
                            } else if(((((twist_ss1 > 93.6) && (twist_ss1 < 103.4)) || ((radius1 > 2.2) && (radius1 < 2.4))) && ((height1 > 0.9) && (height1 < 2.1))) && (twist1 > 180)) {
                                helicalCharacteristics[0] = "a";
                            } else if((((radius1 > 1.1) && (radius1 < 1.9)) && ((twist1 > 223.6) && (twist1 < 253))) && ((height1 > 2.8) && (height1 <= 3.2))) {
                                helicalCharacteristics[0] = "P";
                            }

                            if(((((twist1 > 103.4) && (twist1 < 114.9)) || ((radius1 > 2) && (radius1 <= 2.2)))&&((height1 > 0.9) && (height1 < 2.1))) && (twist1 < 180)) {
                                helicalCharacteristics[1] = "G";
                            } else if(((((twist_ss1 > 103.4) && (twist_ss1 < 114.9)) || ((radius1 > 2) && (radius1 <= 2.2))) && ((height1 > 0.9) && (height1 < 2.1))) && (twist1 > 180)) {
                                helicalCharacteristics[1] = "g";
                            }

                            if(((((twist1 > 77.9) && (twist1 < 93.6)) || ((radius1 > 2.4) && (radius1 <= 2.7))) && ((height1 > 0.9) && (height1 < 2.1))) && (twist1 < 180)) {
                                helicalCharacteristics[2] = "I";
                            } else if(((((twist_ss1 > 77.9) && (twist_ss1 < 93.6)) || ((radius1 > 2.4) && (radius1 <= 2.7))) && ((height1 > 0.9) && (height1 < 2.1))) && (twist1 > 180)) {
                                helicalCharacteristics[2] = "i";
                            }
                        }

                        secondaryStructure1.setCharacteristics(helicalCharacteristics);
                        secondaryStructure1.setStretchId(stretchId);
                    }
                    aminoAcidsInCurrentStretch.stream()
                            .map(this::getSecondaryStructure)
                            .forEach(secondaryStructure -> secondaryStructure.setContinuous(true));
                    stretches.add(new ArrayList<>(aminoAcidsInCurrentStretch));
                    stretchId++;
                }
                aminoAcidsInCurrentStretch.clear();
            }
        }
        return stretches;
    }

    private boolean isContinuousStretch(List<AminoAcid> aminoAcids, int residueIndex) {
        ASSPSecondaryStructure secondaryStructure1 = getSecondaryStructure(aminoAcids.get(residueIndex - 1));
        ASSPSecondaryStructure secondaryStructure2 = getSecondaryStructure(aminoAcids.get(residueIndex));
        return Math.abs(secondaryStructure1.getT() - secondaryStructure2.getT()) <= 35 &&
                Math.abs(secondaryStructure1.getH() - secondaryStructure2.getH()) <= 1.1 &&
                Math.abs(secondaryStructure1.getVtor() - secondaryStructure2.getVtor()) <= 50;
    }

    /**
     * Removes all initial and terminal amino acids whose assignment is 'U U U'.
     * @param stretches the stretches to process
     */
    private void dropInitialAndTerminalCoil(List<List<AminoAcid>> stretches) {
        stretches.forEach(stretch -> {
            List<AminoAcid> aminoAcidsToBeRemoved = new ArrayList<>();
            for(int index = stretch.size() - 1; index >= 0; index--) {
                AminoAcid aminoAcid = stretch.get(index);
                if(getSecondaryStructure(aminoAcid).isUnassigned()) {
                    aminoAcidsToBeRemoved.add(aminoAcid);
                } else {
                    break;
                }
            }
            for(int index = 0; index < stretch.size(); index++) {
                AminoAcid aminoAcid = stretch.get(index);
                if(getSecondaryStructure(stretch.get(index)).isUnassigned()) {
                    aminoAcidsToBeRemoved.add(aminoAcid);
                } else {
                    break;
                }
            }
            stretch.removeAll(aminoAcidsToBeRemoved);
        });
    }

    private void assignFinalCharacteristic(List<List<AminoAcid>> stretches) {
        for (List<AminoAcid> stretch : stretches) {
            // do nothing for stretches less than 2 amino acids
            if(stretch.size() <= 1) {
                continue;
            }

            boolean pp2 = false;
            for(AminoAcid aminoAcid : stretch) {
                ASSPSecondaryStructure secondaryStructure = getSecondaryStructure(aminoAcid);
                String alpha = secondaryStructure.getAlpha();
                String three = secondaryStructure.getThree();
                String pi = secondaryStructure.getPi();
                double radius = secondaryStructure.getR();

                String finalCharacteristics = "U";

                if("A".equalsIgnoreCase(alpha) && "U".equals(three) && "U".equals(pi)) {
                    finalCharacteristics = "A";
                } else if("P".equals(alpha) && "U".equals(three) && "U".equals(pi)) {
                    finalCharacteristics= "P";
                    pp2 = true;
                } else if("A".equalsIgnoreCase(alpha) && "G".equalsIgnoreCase(three) && "U".equals(pi)) {
                    if(radius < 2.2) {
                        finalCharacteristics = "G";
                    } else {
                        finalCharacteristics = "A";
                    }
                } else if("A".equalsIgnoreCase(alpha) && "U".equals(three) && "I".equalsIgnoreCase(pi)) {
                    if(radius > 2.4) {
                        finalCharacteristics = "I";
                    } else {
                        finalCharacteristics = "A";
                    }
                } else if("U".equals(alpha) && "U".equals(three) && "I".equalsIgnoreCase(pi)) {
                    finalCharacteristics = "I";
                } else if("U".equals(alpha) && "G".equalsIgnoreCase(three) && "U".equals(pi)) {
                    finalCharacteristics = "G";
                } else if("A".equalsIgnoreCase(alpha) && "G".equalsIgnoreCase(three) && "I".equalsIgnoreCase(pi)) {
                    if(radius <= 2.2) {
                        finalCharacteristics = "G";
                    } else if(radius > 2.4) {
                        finalCharacteristics = "I";
                    } else {
                        finalCharacteristics = "A";
                    }
                } else if("U".equals(alpha) && "U".equals(three) && "U".equals(pi)) {
                    double twist = secondaryStructure.getT();
                    double mod_t = 360 - twist;
                    double height = secondaryStructure.getH();
                    if((((twist > 104.25) && (twist <= 120)) || ((mod_t > 104.25) && (twist > 180))) && (!pp2)) {
                        finalCharacteristics = "G";
                    } else if(((twist < 92.17) || ((mod_t < 92.17) && (twist  > 180))) && (!pp2)) {
                        finalCharacteristics = "I";
                    } else if((((radius > 1.1) && (radius < 1.9)) || ((twist > 220) && (twist < 260))) && ((height >= 2.7) && (height <= 3.2)) && ((twist > 180))) {
                        finalCharacteristics = "P";
                    } else {
                        if(!pp2 && twist < 120) {
                            finalCharacteristics = "A";
                        } else {
                            finalCharacteristics = "U";
                        }
                    }
                }
                secondaryStructure.setFinalCharacteristic(finalCharacteristics);
            }

            int length = stretch.size();
            String finala = stretch.stream()
                    .map(this::getSecondaryStructure)
                    .map(ASSPSecondaryStructure::getFinalCharacteristic)
                    .collect(Collectors.joining(" "));
            if(length < 3 && (finala.contains("P P P") || ((!finala.contains("P P") && finala.contains("P"))))) {
                finala = "";
            }
            if(length < 2) {
                finala = "";
            }
            if(length == 2 && !finala.contains("G")) {
                if("A G".equals(finala) || "A A".equals(finala) || "G A".equals(finala)) {
                    ASSPSecondaryStructure secondaryStructure1 = getSecondaryStructure(stretch.get(0));
                    ASSPSecondaryStructure secondaryStructure2 = getSecondaryStructure(stretch.get(1));
                    double avg_twist = (secondaryStructure1.getT() + secondaryStructure2.getT()) / 2;
                    double avg_rad = (secondaryStructure1.getR() + secondaryStructure2.getR()) / 2;
                    double mod_t = 360 - avg_twist;
                    if((avg_rad <= 2.3) && (((avg_twist > 102) && (avg_twist <= 120)) || ((mod_t > 102) && (avg_twist > 180)))) {
                        finala = "G G";
                    } else {
                        finala = "";
                    }
                } else {
                    finala = "";
                }
            }
            if(length == 3) {
                if("A G G".equals(finala) || "G A G".equals(finala)) {
                    finala = "G G G";
                }
                if("A G A".equals(finala) || "G A A".equals(finala)) {
                    finala = "A A A";
                }
            }
            if(!"".equals(finala)) {
                finala = streamlineAssignmentString(stretch, finala);
                deriveSecondaryStructureElements(stretch, finala);
            }
        }
    }

    /**
     * ;>
     * @param stretch all residues of this stretch
     * @param finala the string describing the assignments to make
     */
    private String streamlineAssignmentString(List<AminoAcid> stretch, String finala) {
        if(finala.startsWith("A I I I") || finala.startsWith("G I I I") || finala.startsWith("I A I I") || finala.startsWith("A I I A")) {
            finala = "I I I I" + finala.substring(7);
        }
        if(finala.startsWith("A I A A") || finala.startsWith("A G A A") || finala.startsWith("I I A A")) {
            finala = "A A A A" + finala.substring(7);
        }
        if(finala.startsWith("A A I I I") || finala.startsWith("I A A I I")) {
            finala = "I I I I I" + finala.substring(9);
        }
        if(finala.endsWith("I I I A") || finala.endsWith("I I I G") || finala.endsWith("I I A I")) {
            finala = finala.substring(0, finala.length() - 7) + "I I I I";
        }
        if(finala.endsWith("A A I I")) {
            finala = finala.substring(0, finala.length() - 7) + "A A A A";
        }
        if(finala.startsWith("A G G G") || finala.startsWith("G A G G") || finala.startsWith("G G A G") || finala.startsWith("A I G G")) {
            finala = "G G G G" + finala.substring(7);
        }
        if(finala.startsWith("A G G") || finala.startsWith("G A G")) {
            finala = "G G G" + finala.substring(5);
        }
        if(finala.startsWith("I A A") || finala.startsWith("G A A")) {
            finala = "A A A" + finala.substring(5);
        }
        if(finala.endsWith("G G A")) {
            finala = finala.substring(0, finala.length() - 5) + "G G G";
        }
        if(finala.endsWith("A A I") || finala.endsWith("A A G") || finala.endsWith("A G A") || finala.endsWith("A I A")) {
            finala = finala.substring(0, finala.length() - 5) + "A A A";
        }
        if(finala.endsWith("G G G A") || finala.endsWith("G G A G")) {
            finala = finala.substring(0, finala.length() - 7) + "G G G G";
        }
        if(finala.endsWith("G A A G")) {
            finala = finala.substring(0, finala.length() - 7) + "A A A A";
        }
        finala = finala.replace("G G A G G", "G G G G G");
        finala = finala.replace("A A G A A", "A A A A A");
        finala = finala.replace("A A I A A", "A A A A A");
        finala = finala.replace("G G A A G", "G G G G G");
        finala = finala.replace("G A G A G", "G G G G G");
        finala = finala.replace("G A I", "A A A");
        finala = finala.replace("I A G", "A A A");
        finala = finala.replace("I G A", "A A A");
        finala = finala.replace("G I A", "A A A");
        finala = finala.replace("A I A", "A A A");
        finala = finala.replace("A G I", "A A A");
        finala = finala.replace("A I G", "A A A");
        finala = finala.replace("I I A I I", "I I I I I");
        finala = finala.replace("I I A I", "I I I I");
        finala = finala.replace("I I G I", "I I I I");
        finala = finala.replace("A A I U", "A A A U");
        finala = finala.replace("A A G U", "A A A U");
        if(finala.endsWith("A A A I I") || finala.endsWith("A A A G A")) {
            finala = finala.substring(0, finala.length() - 9) + "A A A A A";
        }
        if(finala.endsWith("I I I A A")) {
            finala = finala.substring(0, finala.length() - 9) + "I I I I I";
        }
        return finala;
    }

    private void deriveSecondaryStructureElements(List<AminoAcid> stretch, String finala) {
        String[] split = finala.split(" ");

        boolean inSecondaryStructure = false;
        int startIndex = -1;
        String type = "";

        for(int i = 0; i < split.length; i++) {
            String assignment = split[i];
            if(i < split.length - 1 && assignment.equals(split[i + 1])) {
                inSecondaryStructure = true;
                startIndex = i;
                type = assignment;
            } else {
                if(inSecondaryStructure) {
                    ASSPSecondaryStructure start = getSecondaryStructure(stretch.get(startIndex));
                    ASSPSecondaryStructure end = getSecondaryStructure(stretch.get(i));

                    if(matches(start, "AGI") || (start.isUnassigned() && i != 0 && start.getT() < 180)) {
                        type = type;
                    } else if(matches(start, "agi") || (start.isUnassigned() && i != 0 && start.getT() > 180)) {
                        type = type.toLowerCase();
                    } /*else if(matches(start, "P")) {
                        type = type;
                    }*/
                    if(matches(start, "AGIPagi") || ("AGIagi".contains(type) && i != 0)) {
                        SecondaryStructureElement secondaryStructureElement = matchToSecondaryStructureElement(type);
                        stretch.stream()
                                .map(this::getSecondaryStructure)
                                .forEach(secondaryStructure -> secondaryStructure.setSecondaryStructure(secondaryStructureElement));
                    }

                    inSecondaryStructure = false;
                }
            }
        }
    }

    private SecondaryStructureElement matchToSecondaryStructureElement(String type) {
        switch (type) {
            case "A":case "a":
                return SecondaryStructureElement.ALPHA_HELIX;
            case "P":
                return SecondaryStructureElement.POLYPROLINE_HELIX;
            case "G":case "g":
                return SecondaryStructureElement.THREE_TEN_HELIX;
            case "I":case "i":
                return SecondaryStructureElement.PI_HELIX;
            default:
                return SecondaryStructureElement.COIL;
        }
    }

    private boolean matches(ASSPSecondaryStructure secondaryStructure, String string) {
        for(int i = 0; i < string.length() - 1; i++) {
            String substring = string.substring(i, i + 1);
            if(secondaryStructure.getAlpha().equals(substring) || secondaryStructure.getThree().equals(substring) || secondaryStructure.getPi().equals(substring)) {
                return true;
            }
        }
        return false;
    }

    private ASSPSecondaryStructure getSecondaryStructure(AminoAcid aminoAcid) {
        return aminoAcid.getFeatureContainer().getFeature(ASSPSecondaryStructure.class);
    }
}
