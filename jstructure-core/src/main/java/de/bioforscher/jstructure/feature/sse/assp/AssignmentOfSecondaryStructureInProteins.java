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
        assignHelicalCharacteristics(aminoAcids);
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
                double t = Math.acos(rawTwist) * 57.2958;
                //TODO theta manipulations?

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
                // 57.2958 = 180 / Math.PI, used to convert to degree
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

    private void assignHelicalCharacteristics(List<AminoAcid> aminoAcids) {
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
//                    aminoAcidsInCurrentStretch.add(aminoAcids.get(aminoAcids.indexOf(aminoAcidsInCurrentStretch.get(aminoAcidsInCurrentStretch.size() - 1)) + 1));
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
                    }
                    aminoAcidsInCurrentStretch.stream()
                            .map(this::getSecondaryStructure)
                            .forEach(secondaryStructure -> secondaryStructure.setContinuous(true));
                }
                aminoAcidsInCurrentStretch.clear();
            }
        }

//        // for each fragment of 4 consecutive amino acids: evaluate their helical parameters to derive bend angle
//        for(int residueIndex = 0; residueIndex < aminoAcids.size() - 2; residueIndex++) {
//            AminoAcid aa1 = aminoAcids.get(residueIndex);
//            ASSPSecondaryStructure secondaryStructure1 = getSecondaryStructure(aa1);
//            AminoAcid aa2 = aminoAcids.get(residueIndex + 1);
//            ASSPSecondaryStructure secondaryStructure2 = getSecondaryStructure(aa2);
//            AminoAcid aa3 = aminoAcids.get(residueIndex + 2);
//            ASSPSecondaryStructure secondaryStructure3 = getSecondaryStructure(aa3);
//
//            double twist1 = secondaryStructure1.getT();
//            double twist2 = secondaryStructure2.getT();
//            double twist3 = secondaryStructure3.getT();
//            double height1 = secondaryStructure1.getH();
//            double height2 = secondaryStructure2.getH();
//            double height3 = secondaryStructure3.getH();
//            double radius1 = secondaryStructure1.getR();
//            double radius2 = secondaryStructure2.getR();
//            double radius3 = secondaryStructure3.getR();
//            double sumt = (twist1 + twist2) / 2;
//            double sumh = (height1 + height2) / 2;
//            double sumrad = (radius1 + radius2) / 2;
//            double sumt1 = (twist1 + twist2 + twist3) / 3;
//            double sumh1 = (height1 + height2 + height3) / 3;
//            double sumrad1 = (radius1 + radius2 + radius3) / 3;
//            double twist_ss1 = 360 - twist1;
//            double vtor_ss1 = 360 - secondaryStructure1.getVtor();
//
//            String a1 = "U";
//            String a2 = "U";
//            String a3 = "U";
//
//            double sumt_l = 360 - sumt;
//            double sumt1_l = 360 - sumt1;
//            if(((((sumt > 93.6) && (sumt < 103.5)) || ((twist1 > 93.6) && (twist1 < 103.5)) || ((sumrad > 2.2) && (sumrad < 2.4)) || ((radius1 > 2.2) && (radius1 < 2.4))) && (((sumh > 0.9) && (sumh < 2.1)) || ((height1 > 0.9) && (height1 < 2.1)))) && (twist1 < 180)) {
//                a1 = "A";
//            } else if(((((sumt_l > 93.6) && (sumt_l < 103.5)) || ((twist_ss1 > 93.6) && (twist_ss1 < 103.5)) || ((sumrad > 2.2) && (sumrad < 2.4)) || ((radius1 > 2.2) && (radius1 < 2.4))) && (((sumh > 0.9) && (sumh < 2.1 )) || ((height1 > 0.9) && (height1 < 2.1)))) && (twist1 > 180)) {
//                a1 = "a";
//            } else if((((radius1 > 1.1) && (radius1 < 1.9)) && ((twist1 > 223.6) && (twist1 < 253))) && ((height1 > 2.8) && (height1 <= 3.2))) {
//                a1 = "P";
//            }
//
//            if(((((twist1 > 103.4) && (twist1 < 114.9)) || ((radius1 > 2) && (radius1 <= 2.2))) && ((height1 > 0.9) && (height1 < 2.1))) && (twist1 < 180)) {
//                a2 = "G";
//            } else if(((((twist_ss1 > 103.4) && (twist_ss1 < 114.9)) || ((radius1 > 2) && (radius1 <= 2.2))) && ((height1 > 0.9) && (height1 < 2.1))) && (twist1 > 180)) {
//                a2 = "g";
//            }
//
//            if(((((sumt1 > 77.9) && (sumt1 < 93.6)) || ((sumrad1 > 2.4) && (sumrad1 <= 2.7)) || ((twist1 > 77.9) && (twist1 < 93.6)) || ((radius1 > 2.4) && (radius1 <= 2.7))) && (((sumh1 > 0.9) && (sumh1 < 2.1)) || ((height1 > 0.9) && (height1 < 2.1)))) && (twist1 < 180)) {
//                a3 = "I";
//            } else if(((((sumt1_l > 77.9) && (sumt1_l < 93.6)) || ((sumrad1 > 2.4) && (sumrad1 <= 2.7)) || ((twist_ss1 > 77.9) && (twist_ss1 < 93.6)) || ((radius1 > 2.4) && (radius1 <= 2.7))) && (((sumh1 > 0.9) && (sumh1 < 2.1)) || ((height1 > 0.9) && (height1 < 2.1)))) && (twist1 > 180)) {
//                a3 = "i";
//            }
//
//            secondaryStructure1.setCharacteristics(a1, a2, a3);
//
//            if(Math.abs(twist1 - twist2) <= 35 &&
//                    Math.abs(height1 - height2) <= 1.1 &&
//                    Math.abs(secondaryStructure1.getVtor() - secondaryStructure2.getVtor()) <= 50) {
//                secondaryStructure1.setContinuous(true);
//                secondaryStructure1.setStretchId(stretchId);
//                lastAminoAcidWasContinuous = true;
//            } else {
//                if(lastAminoAcidWasContinuous) {
//                    secondaryStructure1.setContinuous(true);
//                    secondaryStructure1.setStretchId(stretchId);
//                    stretchId++;
//                }
//                lastAminoAcidWasContinuous = false;
//            }
//        }
    }

    private boolean isContinuousStretch(List<AminoAcid> aminoAcids, int residueIndex) {
        ASSPSecondaryStructure secondaryStructure1 = getSecondaryStructure(aminoAcids.get(residueIndex - 1));
        ASSPSecondaryStructure secondaryStructure2 = getSecondaryStructure(aminoAcids.get(residueIndex));
        return Math.abs(secondaryStructure1.getT() - secondaryStructure2.getT()) <= 35 &&
                Math.abs(secondaryStructure1.getH() - secondaryStructure2.getH()) <= 1.1 &&
                Math.abs(secondaryStructure1.getVtor() - secondaryStructure2.getVtor()) <= 50;
    }

    private void assignSecondaryStructureElements(List<AminoAcid> aminoAcids) {
        // for each amino acid: set actually present secondary structure
        for(int residueIndex = 0; residueIndex < aminoAcids.size(); residueIndex++) {
            AminoAcid aminoAcid = aminoAcids.get(residueIndex);
            ASSPSecondaryStructure secondaryStructure = getSecondaryStructure(aminoAcid);

        }
    }

    private ASSPSecondaryStructure getSecondaryStructure(AminoAcid aminoAcid) {
        return aminoAcid.getFeatureContainer().getFeature(ASSPSecondaryStructure.class);
    }
}
