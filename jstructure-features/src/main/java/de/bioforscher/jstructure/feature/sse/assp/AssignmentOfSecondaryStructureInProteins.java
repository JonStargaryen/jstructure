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
import java.util.NoSuchElementException;
import java.util.stream.Collectors;

import static de.bioforscher.jstructure.feature.sse.SecondaryStructureElement.*;

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
        List<RawSecondaryStructure> rawSecondaryStructures = assignFinalCharacteristic(aminoAcids, stretches);

//        List<AminoAcid> filteredAminoAcids = dropInitialAndTerminalCoil(aminoAcids);
//        List<RawSecondaryStructure> rawSecondaryStructures = assignFinalCharacteristic(filteredAminoAcids);

        assignSecondaryStructure(aminoAcids, rawSecondaryStructures);
    }

    /**
     * Assigns the neutral state to an amino acid, i.e. marking its secondary structure as coil.
     * @param aminoAcid the residue to process
     */
    private void assignNeutralState(AminoAcid aminoAcid) {
        aminoAcid.getFeatureContainer().addFeature(new ASSPSecondaryStructure(this, COIL));
    }

    private void assignHelicalParameters(List<AminoAcid> aminoAcids) {
        // for each fragment of 4 consecutive amino acids: evaluate their geometry to derive secondary structure
        for(int residueIndex = 0; residueIndex < aminoAcids.size() - 3; residueIndex++) {
            AminoAcid aa1 = aminoAcids.get(residueIndex);
            AminoAcid aa2 = aminoAcids.get(residueIndex + 1);
            AminoAcid aa3 = aminoAcids.get(residueIndex + 2);
            AminoAcid aa4 = aminoAcids.get(residueIndex + 3);

            if(aa1.getResidueIdentifier().getResidueNumber() + 1 != aa2.getResidueIdentifier().getResidueNumber() ||
                    aa2.getResidueIdentifier().getResidueNumber() + 1 != aa3.getResidueIdentifier().getResidueNumber() ||
                    aa3.getResidueIdentifier().getResidueNumber() + 1 != aa4.getResidueIdentifier().getResidueNumber()) {
                // happens e.g. for residues not solved by xray
                logger.debug("encountered non-consecutive amino acids {} {} {} {}",
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
                logger.debug("missing alpha carbons on fragment {} {} {} {}",
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
                    logger.debug("stretch from {} to {}",
                            aminoAcidsInCurrentStretch.get(1).getResidueIdentifier().getResidueNumber(),
                            aminoAcidsInCurrentStretch.get(aminoAcidsInCurrentStretch.size() - 1).getResidueIdentifier().getResidueNumber() + 2);
                    // if so, determine secondary structure
                    for(int i = 0; i < aminoAcidsInCurrentStretch.size(); i++) {
                        ASSPSecondaryStructure secondaryStructure1 = getSecondaryStructure(aminoAcids.get(aminoAcids.indexOf(aminoAcidsInCurrentStretch.get(i))));
                        ASSPSecondaryStructure secondaryStructure2 = getSecondaryStructure(aminoAcids.get(aminoAcids.indexOf(aminoAcidsInCurrentStretch.get(i)) + 1));
                        double twist1 = secondaryStructure1.getTwist();
                        double twist2 = secondaryStructure2.getTwist();
                        double height1 = secondaryStructure1.getHeight();
                        double height2 = secondaryStructure2.getHeight();
                        double radius1 = secondaryStructure1.getRadius();
                        double radius2 = secondaryStructure2.getRadius();
                        double sumt = (twist1 + twist2) / 2;
                        double sumh = (height1 + height2) / 2;
                        double sumrad = (radius1 + radius2) / 2;
                        double twist_ss1 = 360 - twist1;

                        String[] helicalCharacteristics = { "U", "U", "U" };

                        if(i < aminoAcidsInCurrentStretch.size() - 2) {
                            ASSPSecondaryStructure secondaryStructure3 = getSecondaryStructure(aminoAcidsInCurrentStretch.get(i + 2));
                            double twist3 = secondaryStructure3.getTwist();
                            double sumt1 = (twist1 + twist2 + twist3) / 3;
                            double sumh1 = (height1 + height2 + secondaryStructure3.getHeight()) / 3;
                            double sumrad1 = (radius1 + radius2 + secondaryStructure3.getRadius()) / 3;
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
                            //TODO deviation for ILE-30: expected 'A U U', actually 'A U I'
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

                        //TODO remove - temporary fix for other bug fixes downstream - did not cause problems
//                        int resNum = aminoAcids.get(aminoAcids.indexOf(aminoAcidsInCurrentStretch.get(i))).getResidueIdentifier().getResidueIdentifier();
//                        if(resNum == 30 || resNum == 33) {
//                            secondaryStructure1.setCharacteristics(new String[] { "A", "U", "U" });
//                        }
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
        return Math.abs(secondaryStructure1.getTwist() - secondaryStructure2.getTwist()) <= 35 &&
                Math.abs(secondaryStructure1.getHeight() - secondaryStructure2.getHeight()) <= 1.1 &&
                Math.abs(secondaryStructure1.getVirtualTorsionAngle() - secondaryStructure2.getVirtualTorsionAngle()) <= 50;
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

//    private List<AminoAcid> dropInitialAndTerminalCoil(List<AminoAcid> aminoAcids) {
//        return aminoAcids.stream()
//                .filter(aminoAcid -> {
//                    ASSPSecondaryStructure secondaryStructure = getSecondaryStructure(aminoAcid);
//                    try {
//                        return !secondaryStructure.isUnassigned();
//                    } catch (NullPointerException e) {
//                        return false;
//                    }
//                })
//                .collect(Collectors.toList());
//    }

    private List<RawSecondaryStructure> assignFinalCharacteristic(List<AminoAcid> aminoAcids, List<List<AminoAcid>> stretches) {
        List<RawSecondaryStructure> rawSecondaryStructures = new ArrayList<>();
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
                double radius = secondaryStructure.getRadius();

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
                    double twist = secondaryStructure.getTwist();
                    double mod_t = 360 - twist;
                    double height = secondaryStructure.getHeight();
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
                    double avg_twist = (secondaryStructure1.getTwist() + secondaryStructure2.getTwist()) / 2;
                    double avg_rad = (secondaryStructure1.getRadius() + secondaryStructure2.getRadius()) / 2;
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
                finala = streamlineAssignmentString(finala);
                List<RawSecondaryStructure> rss = deriveRawSecondaryStructureElements(stretch, finala);
                rawSecondaryStructures.addAll(reorganize(aminoAcids, rss));
            }
        }

        return rawSecondaryStructures;
    }

//    private List<RawSecondaryStructure> assignFinalCharacteristic(List<AminoAcid> aminoAcids) {
//        List<RawSecondaryStructure> rawSecondaryStructures = new ArrayList<>();
//        boolean pp2 = false;
//        for(AminoAcid aminoAcid : aminoAcids) {
//            ASSPSecondaryStructure secondaryStructure = getSecondaryStructure(aminoAcid);
//            String alpha = secondaryStructure.getAlpha();
//            String three = secondaryStructure.getThree();
//            String pi = secondaryStructure.getPi();
//            double radius = secondaryStructure.getRadius();
//
//            String finalCharacteristics = "U";
//
//            if("A".equalsIgnoreCase(alpha) && "U".equals(three) && "U".equals(pi)) {
//                finalCharacteristics = "A";
//            } else if("P".equals(alpha) && "U".equals(three) && "U".equals(pi)) {
//                finalCharacteristics= "P";
//                pp2 = true;
//            } else if("A".equalsIgnoreCase(alpha) && "G".equalsIgnoreCase(three) && "U".equals(pi)) {
//                if(radius < 2.2) {
//                    finalCharacteristics = "G";
//                } else {
//                    finalCharacteristics = "A";
//                }
//            } else if("A".equalsIgnoreCase(alpha) && "U".equals(three) && "I".equalsIgnoreCase(pi)) {
//                if(radius > 2.4) {
//                    finalCharacteristics = "I";
//                } else {
//                    finalCharacteristics = "A";
//                }
//            } else if("U".equals(alpha) && "U".equals(three) && "I".equalsIgnoreCase(pi)) {
//                finalCharacteristics = "I";
//            } else if("U".equals(alpha) && "G".equalsIgnoreCase(three) && "U".equals(pi)) {
//                finalCharacteristics = "G";
//            } else if("A".equalsIgnoreCase(alpha) && "G".equalsIgnoreCase(three) && "I".equalsIgnoreCase(pi)) {
//                if(radius <= 2.2) {
//                    finalCharacteristics = "G";
//                } else if(radius > 2.4) {
//                    finalCharacteristics = "I";
//                } else {
//                    finalCharacteristics = "A";
//                }
//            } else if("U".equals(alpha) && "U".equals(three) && "U".equals(pi)) {
//                double twist = secondaryStructure.getTwist();
//                double mod_t = 360 - twist;
//                double height = secondaryStructure.getHeight();
//                if((((twist > 104.25) && (twist <= 120)) || ((mod_t > 104.25) && (twist > 180))) && (!pp2)) {
//                    finalCharacteristics = "G";
//                } else if(((twist < 92.17) || ((mod_t < 92.17) && (twist  > 180))) && (!pp2)) {
//                    finalCharacteristics = "I";
//                } else if((((radius > 1.1) && (radius < 1.9)) || ((twist > 220) && (twist < 260))) && ((height >= 2.7) && (height <= 3.2)) && ((twist > 180))) {
//                    finalCharacteristics = "P";
//                } else {
//                    if(!pp2 && twist < 120) {
//                        finalCharacteristics = "A";
//                    } else {
//                        finalCharacteristics = "U";
//                    }
//                }
//            }
//            secondaryStructure.setFinalCharacteristic(finalCharacteristics);
//        }
//
//        int length = aminoAcids.size();
//        String finala = aminoAcids.stream()
//                .map(this::getSecondaryStructure)
//                .map(ASSPSecondaryStructure::getFinalCharacteristic)
//                .collect(Collectors.joining(" "));
//        if(length < 3 && (finala.contains("P P P") || ((!finala.contains("P P") && finala.contains("P"))))) {
//            finala = "";
//        }
//        if(length < 2) {
//            finala = "";
//        }
//        if(length == 2 && !finala.contains("G")) {
//            if("A G".equals(finala) || "A A".equals(finala) || "G A".equals(finala)) {
//                ASSPSecondaryStructure secondaryStructure1 = getSecondaryStructure(aminoAcids.get(0));
//                ASSPSecondaryStructure secondaryStructure2 = getSecondaryStructure(aminoAcids.get(1));
//                double avg_twist = (secondaryStructure1.getTwist() + secondaryStructure2.getTwist()) / 2;
//                double avg_rad = (secondaryStructure1.getRadius() + secondaryStructure2.getRadius()) / 2;
//                double mod_t = 360 - avg_twist;
//                if((avg_rad <= 2.3) && (((avg_twist > 102) && (avg_twist <= 120)) || ((mod_t > 102) && (avg_twist > 180)))) {
//                    finala = "G G";
//                } else {
//                    finala = "";
//                }
//            } else {
//                finala = "";
//            }
//        }
//        if(length == 3) {
//            if("A G G".equals(finala) || "G A G".equals(finala)) {
//                finala = "G G G";
//            }
//            if("A G A".equals(finala) || "G A A".equals(finala)) {
//                finala = "A A A";
//            }
//        }
//        if(!"".equals(finala)) {
//            finala = streamlineAssignmentString(finala);
//            List<RawSecondaryStructure> rss = deriveRawSecondaryStructureElements(aminoAcids, finala);
//            rawSecondaryStructures.addAll(reorganize(aminoAcids, rss));
//        }
//
//        return rawSecondaryStructures;
//    }

    /**
     * ;>
     * @param finala the string describing the assignments to make
     */
    private String streamlineAssignmentString(String finala) {
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

    private List<RawSecondaryStructure> deriveRawSecondaryStructureElements(List<AminoAcid> aminoAcids, String finala) {
        List<RawSecondaryStructure> secondaryStructures = new ArrayList<>();
        String[] split = finala.split(" ");

        boolean inSecondaryStructure = false;
        int startIndex = -1;
        String type = "";

        for(int i = 0; i < split.length; i++) {
            String assignment = split[i];
            if("U".equals(assignment)) {
                continue;
            }
            if(i < split.length - 1 && assignment.equals(split[i + 1])) {
                inSecondaryStructure = true;
                if(startIndex == -1) {
                    startIndex = i;
                }
                type = assignment;
            } else {
                if(inSecondaryStructure) {
                    ASSPSecondaryStructure start = getSecondaryStructure(aminoAcids.get(startIndex));

                    if(matches(start, "agi") || (start.isUnassigned() && i != 0 && start.getTwist() > 180)) {
                        type = type.toLowerCase();
                    }
                    if(matches(start, "AGIPagi") || ("AGIagi".contains(type) && i != 0)) {
                        SecondaryStructureElement secondaryStructureElement = matchToSecondaryStructureElement(type);
                        secondaryStructures.add(new RawSecondaryStructure(secondaryStructureElement,
                                type,
                                aminoAcids.get(startIndex).getResidueIdentifier().getResidueNumber() + 1,
                                aminoAcids.get(i).getResidueIdentifier().getResidueNumber() + 2));
                    }

                    startIndex = -1;
                    inSecondaryStructure = false;
                }
            }
        }
        return secondaryStructures;
    }

    class RawSecondaryStructure {
        SecondaryStructureElement secondaryStructureElement;
        String type;
        int start;
        int end;
        int length;

        RawSecondaryStructure(SecondaryStructureElement secondaryStructureElement,
                              String type,
                              int start,
                              int end) {
            this.secondaryStructureElement = secondaryStructureElement;
            this.type = type;
            this.start = start;
            this.end = end;
            this.length = end - start + 1;
        }

        @Override
        public String toString() {
            return secondaryStructureElement + "\t" + start + "\t" + end + "\t" + length;
        }

        boolean sameNatureAs(RawSecondaryStructure other) {
            return isFreakSecondaryStructure() == other.isFreakSecondaryStructure();
        }

        boolean isFreakSecondaryStructure() {
            return !("A".equals(type) || "G".equals(type) || "I".equals(type) || "S".equals(type));
        }

        boolean isMinimalLength() {
            if(secondaryStructureElement == ALPHA_HELIX) {
                return length == 4;
            } else if(secondaryStructureElement == THREE_TEN_HELIX || secondaryStructureElement == POLYPROLINE_HELIX) {
                return length == 3;
            } else if(secondaryStructureElement == PI_HELIX) {
                return length == 5;
            } else {
                return false;
            }
        }
    }

    private SecondaryStructureElement matchToSecondaryStructureElement(String type) {
        switch (type) {
            case "A":case "a":
                return ALPHA_HELIX;
            case "P":
                return POLYPROLINE_HELIX;
            case "G":case "g":
                return THREE_TEN_HELIX;
            case "I":case "i":
                return PI_HELIX;
            default:
                return COIL;
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

    private List<RawSecondaryStructure> reorganize(List<AminoAcid> aminoAcids, List<RawSecondaryStructure> rawSecondaryStructures) {
        List<RawSecondaryStructure> rawSecondaryStructuresToRemove = new ArrayList<>();

        for(int i = 1; i < rawSecondaryStructures.size() - 1; i++) {
            RawSecondaryStructure rawSecondaryStructurePrevious = rawSecondaryStructures.get(i - 1);
            RawSecondaryStructure rawSecondaryStructure = rawSecondaryStructures.get(i);
            RawSecondaryStructure rawSecondaryStructureNext = rawSecondaryStructures.get(i + 1);

            if((("a".equals(rawSecondaryStructure.type)) ||
                    ("P".equals(rawSecondaryStructure.type)) ||
                    ("i".equals(rawSecondaryStructure.type)) ||
                    ("g".equals(rawSecondaryStructure.type))) &&
                    (((!"a".equals(rawSecondaryStructurePrevious.type)) ||
                    (!"P".equals(rawSecondaryStructurePrevious.type)) ||
                    (!"i".equals(rawSecondaryStructurePrevious.type)) ||
                    (!"g".equals(rawSecondaryStructurePrevious.type))))) {
                if(rawSecondaryStructure.end >= rawSecondaryStructureNext.start) {
                    if(("a".equals(rawSecondaryStructure.type) && rawSecondaryStructure.length < 5) ||
                            ("i".equals(rawSecondaryStructure.type) && rawSecondaryStructure.length < 6) ||
                            ("g".equals(rawSecondaryStructure.type) && rawSecondaryStructure.length < 4) ||
                            ("P".equals(rawSecondaryStructure.type) && rawSecondaryStructure.length < 4)) {
                        rawSecondaryStructuresToRemove.add(rawSecondaryStructure);
                    } else {
//                        logger.info("decreasing length of {} wrt {}", rawSecondaryStructure, rawSecondaryStructureNext);
                        rawSecondaryStructure.end--;
                        rawSecondaryStructure.length--;
                        //TODO stretch ids
                    }
                }
                if(rawSecondaryStructure.start <= rawSecondaryStructurePrevious.end) {
                    if(("a".equals(rawSecondaryStructure.type) && rawSecondaryStructure.length < 5) ||
                            ("i".equals(rawSecondaryStructure.type) && rawSecondaryStructure.length < 6) ||
                            ("g".equals(rawSecondaryStructure.type) && rawSecondaryStructure.length < 4) ||
                            ("P".equals(rawSecondaryStructure.type) && rawSecondaryStructure.length < 4)) {
                        rawSecondaryStructuresToRemove.add(rawSecondaryStructure);
                    } else {
//                        logger.info("decreasing length of {} wrt {}", rawSecondaryStructure, rawSecondaryStructurePrevious);
                        rawSecondaryStructure.start++;
                        rawSecondaryStructure.length--;

                    }
                }
            } else if(rawSecondaryStructure.type.equals(rawSecondaryStructureNext.type) && rawSecondaryStructure.end == rawSecondaryStructureNext.start) {
                AminoAcid aminoAcid = findAminoAcid(aminoAcids, rawSecondaryStructure.end);
                if(getSecondaryStructure(aminoAcid).getBendAngle() < 60) {
//                    logger.info("fusing {} onto {}", rawSecondaryStructureNext, rawSecondaryStructure);
                    rawSecondaryStructure.end = rawSecondaryStructureNext.end;
                    rawSecondaryStructure.length += rawSecondaryStructureNext.length - 1;
                    rawSecondaryStructuresToRemove.add(rawSecondaryStructureNext);
                } else {
                    rawSecondaryStructuresToRemove.addAll(reorganizePair(aminoAcids, rawSecondaryStructure, rawSecondaryStructureNext));
                }
            } else if(!rawSecondaryStructure.type.equals(rawSecondaryStructureNext.type) && rawSecondaryStructure.end == rawSecondaryStructureNext.start) {
                rawSecondaryStructuresToRemove.addAll(reorganizePair(aminoAcids, rawSecondaryStructure, rawSecondaryStructureNext));
            }
        }
        rawSecondaryStructures.removeAll(rawSecondaryStructuresToRemove);
        return rawSecondaryStructures;
    }

    private AminoAcid findAminoAcid(List<AminoAcid> aminoAcids, int residueNumber) {
        for(AminoAcid aminoAcid : aminoAcids) {
            if (aminoAcid.getResidueIdentifier().getResidueNumber() == residueNumber) {
                return aminoAcid;
            }
        }
        throw new NoSuchElementException("did not find amino acid with residue number '" + residueNumber + "'");
    }

    private List<RawSecondaryStructure> reorganizePair(List<AminoAcid> aminoAcids, RawSecondaryStructure rawSecondaryStructure1, RawSecondaryStructure rawSecondaryStructure2) {
        List<RawSecondaryStructure> rawSecondaryStructuresToRemove = new ArrayList<>();

        SecondaryStructureElement secondaryStructureElement1 = rawSecondaryStructure1.secondaryStructureElement;
        SecondaryStructureElement secondaryStructureElement2 = rawSecondaryStructure2.secondaryStructureElement;

        if(((secondaryStructureElement2 == ALPHA_HELIX && rawSecondaryStructure2.length > 4) ||
                (secondaryStructureElement2 == PI_HELIX && rawSecondaryStructure2.length > 5) ||
                (secondaryStructureElement2 == THREE_TEN_HELIX && rawSecondaryStructure2.length > 3)) &&
                ((secondaryStructureElement1 == ALPHA_HELIX && rawSecondaryStructure1.length > 4) ||
                (secondaryStructureElement1 == PI_HELIX && rawSecondaryStructure1.length > 5) ||
                (secondaryStructureElement1 == THREE_TEN_HELIX && rawSecondaryStructure1.length > 3)) &&
                (rawSecondaryStructure1.sameNatureAs(rawSecondaryStructure2))) {
            if (secondaryStructureElement1 == PI_HELIX || secondaryStructureElement1 == THREE_TEN_HELIX || rawSecondaryStructure1.type.equals(rawSecondaryStructure2.type)) {
//                logger.info("[1] decreasing length of {} wrt {}", rawSecondaryStructure2, rawSecondaryStructure1);
                rawSecondaryStructure2.start++;
                rawSecondaryStructure2.length--;
            }
            if(secondaryStructureElement1 == ALPHA_HELIX && (secondaryStructureElement2 == PI_HELIX || secondaryStructureElement2 == THREE_TEN_HELIX)) {
                ASSPSecondaryStructure secondaryStructure1 = getSecondaryStructure(findAminoAcid(aminoAcids, rawSecondaryStructure1.end));
                ASSPSecondaryStructure secondaryStructure2 = getSecondaryStructure(findAminoAcid(aminoAcids, rawSecondaryStructure1.end + 1));
                String alpha1 = secondaryStructure1.getAlpha();
                String three1 = secondaryStructure1.getThree();
                String pi1 = secondaryStructure1.getPi();
                String alpha2 = secondaryStructure2.getAlpha();
                String three2 = secondaryStructure2.getThree();
                String pi2 = secondaryStructure2.getPi();
                if ((("I".equals(alpha1) || "I".equals(three1) || "I".equals(pi1)) && ("I".equals(alpha2) || "I".equals(three2) || "I".equals(pi2)) && secondaryStructureElement2 == PI_HELIX) ||
                        (("G".equals(alpha1) || "G".equals(three1) || "I".equals(pi1)) && ("G".equals(alpha2) || "I".equals(three2) || "I".equals(pi2)) && secondaryStructureElement2 == THREE_TEN_HELIX)) {
//                    logger.info("[2] decreasing length of {} wrt {}", rawSecondaryStructure1, rawSecondaryStructure2);
                    rawSecondaryStructure1.end--;
                    rawSecondaryStructure1.length--;
                } else {
//                    logger.info("[3] decreasing length of {} wrt {}", rawSecondaryStructure2, rawSecondaryStructure2);
                    rawSecondaryStructure2.start++;
                    rawSecondaryStructure2.length--;
                }
            }
        } else if(((secondaryStructureElement2 == ALPHA_HELIX && rawSecondaryStructure2.length > 4) || (secondaryStructureElement2 == PI_HELIX && rawSecondaryStructure2.length > 5) || (secondaryStructureElement2 == THREE_TEN_HELIX && rawSecondaryStructure2.length > 3)) && ((secondaryStructureElement1 == ALPHA_HELIX && rawSecondaryStructure1.length == 4) || (secondaryStructureElement1 == PI_HELIX && rawSecondaryStructure1.length == 5) || (secondaryStructureElement1 == THREE_TEN_HELIX && rawSecondaryStructure1.length == 3))) {
//            logger.info("[4] decreasing length of {} wrt {}", rawSecondaryStructure2, rawSecondaryStructure1);
            rawSecondaryStructure2.start++;
            rawSecondaryStructure2.length--;
        } else if(((secondaryStructureElement2 == ALPHA_HELIX && rawSecondaryStructure2.length == 4) || (secondaryStructureElement2 == PI_HELIX && rawSecondaryStructure2.length == 5) || (secondaryStructureElement2 == THREE_TEN_HELIX && rawSecondaryStructure2.length == 3)) && ((secondaryStructureElement1 == ALPHA_HELIX && rawSecondaryStructure1.length > 4) || (secondaryStructureElement1 == PI_HELIX  && rawSecondaryStructure1.length > 5) || (secondaryStructureElement1 == THREE_TEN_HELIX && rawSecondaryStructure1.length > 3))) {
//            logger.info("[5] decreasing length of {} wrt {}", rawSecondaryStructure1, rawSecondaryStructure2);
            rawSecondaryStructure1.end--;
            rawSecondaryStructure1.length--;
        } else if(((secondaryStructureElement2 == ALPHA_HELIX && rawSecondaryStructure2.length == 4) || (secondaryStructureElement2 == PI_HELIX && rawSecondaryStructure2.length == 5) || (secondaryStructureElement2 == THREE_TEN_HELIX && rawSecondaryStructure2.length == 3)) && ((secondaryStructureElement1 == ALPHA_HELIX && rawSecondaryStructure1.length == 4) || (secondaryStructureElement1 == PI_HELIX && rawSecondaryStructure1.length == 5) || (secondaryStructureElement1 == THREE_TEN_HELIX && rawSecondaryStructure1.length == 3))) {
            if(rawSecondaryStructure1.type.equals(rawSecondaryStructure2.type)) {
//                logger.info("[6] fusing {} onto {}", rawSecondaryStructure2, rawSecondaryStructure1);
                rawSecondaryStructure1.end = rawSecondaryStructure2.end;
                rawSecondaryStructure1.length += rawSecondaryStructure2.length - 1;
                rawSecondaryStructuresToRemove.add(rawSecondaryStructure2);
            } else if(rawSecondaryStructure1.length > rawSecondaryStructure2.length) {
//                logger.info("[7] fusing {} onto {}", rawSecondaryStructure2, rawSecondaryStructure1);
                rawSecondaryStructure1.end = rawSecondaryStructure2.end;
                rawSecondaryStructure1.length += rawSecondaryStructure2.length - 1;
                rawSecondaryStructuresToRemove.add(rawSecondaryStructure2);
            } else if(rawSecondaryStructure2.length > rawSecondaryStructure1.length) {
//                logger.info("[8] fusing {} onto {}", rawSecondaryStructure1, rawSecondaryStructure2);
                rawSecondaryStructure2.length += rawSecondaryStructure1.length - 1;
                rawSecondaryStructure2.start = rawSecondaryStructure1.start;
                rawSecondaryStructuresToRemove.add(rawSecondaryStructure1);
            }
        } else if(((secondaryStructureElement2 == ALPHA_HELIX && rawSecondaryStructure2.length < 4) || (secondaryStructureElement2 == PI_HELIX && rawSecondaryStructure2.length < 5) || (secondaryStructureElement2 == THREE_TEN_HELIX && rawSecondaryStructure2.length < 3)) && ((secondaryStructureElement1 == ALPHA_HELIX && rawSecondaryStructure1.length >= 4) || (secondaryStructureElement1 == PI_HELIX && rawSecondaryStructure1.length >= 5) || (secondaryStructureElement1 == THREE_TEN_HELIX && rawSecondaryStructure1.length >= 3))) {
//            logger.info("[9] fusing {} onto {}", rawSecondaryStructure2, rawSecondaryStructure1);
            rawSecondaryStructure1.end = rawSecondaryStructure2.end;
            rawSecondaryStructure1.length += rawSecondaryStructure2.length - 1;
            rawSecondaryStructuresToRemove.add(rawSecondaryStructure2);
        } else if(((secondaryStructureElement2 == ALPHA_HELIX && rawSecondaryStructure2.length >= 4) || (secondaryStructureElement2 == PI_HELIX && rawSecondaryStructure2.length >= 5) || (secondaryStructureElement2 == THREE_TEN_HELIX && rawSecondaryStructure2.length >= 3)) && ((secondaryStructureElement1 == ALPHA_HELIX && rawSecondaryStructure1.length < 4) || (secondaryStructureElement1 == PI_HELIX && rawSecondaryStructure1.length < 5) || (secondaryStructureElement1 == THREE_TEN_HELIX && rawSecondaryStructure1.length < 3 ))) {
//            logger.info("[10] fusing {} onto {}", rawSecondaryStructure1, rawSecondaryStructure2);
            rawSecondaryStructure2.length += rawSecondaryStructure1.length - 1;
            rawSecondaryStructure2.start = rawSecondaryStructure1.start;
            rawSecondaryStructuresToRemove.add(rawSecondaryStructure1);
        }

        return rawSecondaryStructuresToRemove;
    }

    private void assignSecondaryStructure(List<AminoAcid> aminoAcids, List<RawSecondaryStructure> rawSecondaryStructures) {
        List<RawSecondaryStructure> rawSecondaryStructuresToRemove = new ArrayList<>();
        List<RawSecondaryStructure> alphaHelices = new ArrayList<>();
        boolean pi2alpha = false;

        for (RawSecondaryStructure rawSecondaryStructure : rawSecondaryStructures) {
            int length = rawSecondaryStructure.length;
            String type = rawSecondaryStructure.type;
            if (length < 3) {
                rawSecondaryStructuresToRemove.add(rawSecondaryStructure);
                continue;
            }
            if (rawSecondaryStructure.secondaryStructureElement == ALPHA_HELIX && length > 3) {
                if ("A".equals(type)) {
                    if (!alphaHelices.isEmpty()) {
                        RawSecondaryStructure previousAlphaHelix = alphaHelices.get(alphaHelices.size() - 1);
                        if (pi2alpha && previousAlphaHelix.end == rawSecondaryStructure.start - 1) {
                            rawSecondaryStructure.start = previousAlphaHelix.start;
                            rawSecondaryStructure.length += previousAlphaHelix.length;
                            rawSecondaryStructuresToRemove.add(previousAlphaHelix);
                            alphaHelices.remove(previousAlphaHelix);
                            pi2alpha = false;
                        }
                    }

                    alphaHelices.add(rawSecondaryStructure);
                }
            }
            if (rawSecondaryStructure.secondaryStructureElement == THREE_TEN_HELIX && length > 2) {
                if ("A".equals(type) || "G".equals(type)) {
                    rawSecondaryStructure.type = "G";
                }
                if ("a".equals(type) || "g".equals(type)) {
                    rawSecondaryStructure.type = "g";
                }
            }
            if (rawSecondaryStructure.secondaryStructureElement == PI_HELIX && length > 4) {
                if ("I".equals(type)) {
                    rawSecondaryStructure.secondaryStructureElement = piTwistCheck(aminoAcids, rawSecondaryStructure);
                    if (rawSecondaryStructure.secondaryStructureElement == ALPHA_HELIX) {
                        if (!alphaHelices.isEmpty()) {
                            RawSecondaryStructure previousAlphaHelix = alphaHelices.get(alphaHelices.size() - 1);
                            if (previousAlphaHelix.end == rawSecondaryStructure.start - 1) {
                                pi2alpha = true;
                                rawSecondaryStructuresToRemove.add(previousAlphaHelix);
                                alphaHelices.remove(previousAlphaHelix);
                            }
                        }

                        alphaHelices.add(rawSecondaryStructure);
                    }
                }
            }
        }

        rawSecondaryStructures.removeAll(rawSecondaryStructuresToRemove);

//        naiveCheckForUnassignedRegions(aminoAcids, rawSecondaryStructures);
        naiveCheckForOverlaps(rawSecondaryStructures);

        // finally assign the secondary structure elements - was kinda hard to get here
        rawSecondaryStructures.forEach(rawSecondaryStructure -> {
            for(int i = rawSecondaryStructure.start; i <= rawSecondaryStructure.end; i++) {
                getSecondaryStructure(aminoAcids.get(i)).setSecondaryStructure(rawSecondaryStructure.secondaryStructureElement);
            }
        });
    }

//    private void naiveCheckForUnassignedRegions(List<AminoAcid> aminoAcids, List<RawSecondaryStructure> rawSecondaryStructures) {
//        List<RawSecondaryStructure> rawSecondaryStructuresToRemove = new ArrayList<>();
//        for(RawSecondaryStructure rawSecondaryStructure : rawSecondaryStructures) {
//            for(int i = rawSecondaryStructure.start; i <= rawSecondaryStructure.end; i++) {
//                ASSPSecondaryStructure asspSecondaryStructure = getSecondaryStructure(aminoAcids.get(i));
//                if(!asspSecondaryStructure.isContinuous()) {
//                    System.out.println(aminoAcids.get(i));
//                    System.out.println(asspSecondaryStructure.isContinuous());
//                    System.out.println(asspSecondaryStructure);
//                }
//            }
//        }
//        rawSecondaryStructures.removeAll(rawSecondaryStructuresToRemove);
//    }

    private void naiveCheckForOverlaps(List<RawSecondaryStructure> rawSecondaryStructures) {
        List<RawSecondaryStructure> rawSecondaryStructuresToRemove = new ArrayList<>();
//        for(RawSecondaryStructure rss1 : rawSecondaryStructures) {
//            for(RawSecondaryStructure rss2 : rawSecondaryStructures) {
//                if(rss1 == rss2) {
//                    continue;
//                }
//                if(rss1.end == rss2.start) {
////                    System.out.println(rss2);
//                    if(rss1.type.equals(rss2.type)) {
//                        rss1.end = rss2.end;
//                        rss1.length += rss2.length - 1;
//                        rawSecondaryStructuresToRemove.add(rss2);
//                    } else {
//                        rss2.start++;
//                        rss2.length--;
//                    }
//                } else if(rss1.end + 1 == rss2.start && rss1.type.equals(rss2.type)) {
////                    System.out.println(rss2);
//                    rss1.end = rss2.end;
//                    rss1.length += rss2.length;
//                    rawSecondaryStructuresToRemove.add(rss2);
//                }
//            }
//        }
        for(int i = 0; i < rawSecondaryStructures.size() - 1; i++) {
            RawSecondaryStructure rss1 = rawSecondaryStructures.get(i);
            for(int j = i + 1; j < rawSecondaryStructures.size(); j++) {
                RawSecondaryStructure rss2 = rawSecondaryStructures.get(j);
                if(rss1.end == rss2.start) {
//                    System.out.println(rss2);
                    if(rss1.type.equals(rss2.type)) {
                        rss1.end = rss2.end;
                        rss1.length += rss2.length - 1;
                        rawSecondaryStructuresToRemove.add(rss2);
                    } else {
                        rss2.start++;
                        rss2.length--;
                    }
                } else if(rss1.end + 1 == rss2.start && rss1.type.equals(rss2.type)) {
//                    System.out.println(rss2);
                    rss1.end = rss2.end;
                    rss1.length += rss2.length;
                    rawSecondaryStructuresToRemove.add(rss2);
                }
            }
        }
        rawSecondaryStructures.removeAll(rawSecondaryStructuresToRemove);
        if(!rawSecondaryStructuresToRemove.isEmpty()) {
            naiveCheckForOverlaps(rawSecondaryStructures);
        }
    }

    private SecondaryStructureElement piTwistCheck(List<AminoAcid> aminoAcids, RawSecondaryStructure rawSecondaryStructure) {
        double tw = 0;
        double rad = 0;
        for(AminoAcid aminoAcid : aminoAcids) {
            int resNum = aminoAcid.getResidueIdentifier().getResidueNumber();
            if(resNum >= rawSecondaryStructure.start - 1 && resNum + 3 <= rawSecondaryStructure.end + 1) {
                ASSPSecondaryStructure asspSecondaryStructure = getSecondaryStructure(aminoAcid);
                tw += asspSecondaryStructure.getTwist();
                rad += asspSecondaryStructure.getRadius();
            }
        }

        // normalize by length
        tw /= rawSecondaryStructure.length;
        rad /= rawSecondaryStructure.length;

        if(tw < 93.8 && rad > 2.45) {
            return PI_HELIX;
        } else if(tw > 102 && rad < 2.3) {
            return THREE_TEN_HELIX;
        } else {
            return ALPHA_HELIX;
        }
    }

    private ASSPSecondaryStructure getSecondaryStructure(AminoAcid aminoAcid) {
        return aminoAcid.getFeatureContainer().getFeature(ASSPSecondaryStructure.class);
    }
}
