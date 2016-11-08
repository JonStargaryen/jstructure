package de.bioforscher.jstructure.feature.sse;

import de.bioforscher.jstructure.feature.FeatureProvider;
import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.model.Fragment;
import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.structure.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.stream.Collectors;

/**
 * Annotates the secondary structure element of each getResidue in a
 * {@link Protein}.<br />
 * This is BioJava-code which itself it strongly motivated by Kabsch & Sander's
 * DSSP.<br />
 * original doc:<br />
 * <br />
 *
 * Calculate and assign the secondary structure (SS) to the Groups of a
 * Structure object. This object also stores the result of the calculation.
 * <p>
 * The rules for SS calculation are the ones defined by DSSP: Kabsch,W. and
 * Sander,C. (1983) Biopolymers 22, 2577-2637. Original DSSP article see at:
 * <a href="http://www.cmbi.kun.nl/gv/dssp/dssp.pdf">dssp.pdf</a>. Some parts
 * are also taken from: T.E.Creighton, Proteins - Structure and Molecular
 * Properties, 2nd Edition, Freeman 1994.
 *
 * @author Andreas Prlic
 * @author Aleix Lafita
 * @author Anthony Bradley
 *
 */
public class SecondaryStructureAnnotator implements FeatureProvider {
    final Logger logger = LoggerFactory.getLogger(SecondaryStructureAnnotator.class);
    /**
     * DSSP assigns helices one getResidue shorter at each end, because the
     * residues at (i-1) and (i+n+1) are not assigned helix type although they
     * contain a consistent turn (H-bond). If this parameter is true, the
     * helices will be the length of the original DSSP convention. If it is
     * false, they will be two getResidue longer.
     */
    private static final boolean DSSP_HELICES = true;

    /** min distance between two residues */
    public static final double MINDIST = 0.5;

    /** min distance of two CA atoms if H-bonds are allowed to form */
    public static final double CA_MIN_DIST = 9.0;

    /** max distance CA atoms in peptide bond (backbone discontinuity) */
    public static final double MAX_PEPTIDE_BOND_LENGTH = 2.5;
    /** squared value for faster computations */
    public static final double MAX_PEPTIDE_BOND_LENGTH_SQUARED = MAX_PEPTIDE_BOND_LENGTH * MAX_PEPTIDE_BOND_LENGTH;

    /** Minimal H-bond energy in cal/mol */
    public static final int HBONDLOWENERGY = -9900;

    /** higher limit for H-bond energy */
    public static final double HBONDHIGHENERGY = -500.0;

    /**
     * constant for electrostatic energy
     *
     * <pre>
     *      f  *  q1 *   q2  *  scale
     * Q = -332 * 0.42 * 0.20 * 1000.0
     * </pre>
     *
     * q1 and q2 are partial charges which are placed on the C,O (+q1,-q1) and
     * N,H (-q2,+q2)
     */
    public static final double Q = -27888.0;

    private List<BetaBridge> bridges;
    private List<Ladder> ladders;

    public enum FeatureNames {
        SECONDARY_STRUCTURE_STATES
    }

    private List<Residue> residues;

    @Override
    public void process(Protein protein) {
        residues = protein.residues().collect(Collectors.toList());
        bridges = new ArrayList<>();
        ladders = new ArrayList<>();
        // init mapping
        residues.forEach(r -> r.setFeature(FeatureNames.SECONDARY_STRUCTURE_STATES,
                new SecStrucState(DSSPSecondaryStructureElement.COIL)));

        calculateHAtoms();
        calculateHBonds();
        calculateDihedralAngles();
        calculateTurns();
        buildHelices();
        detectBends();
        detectStrands();

        protein.clearPseudoAtoms();
//        residues.stream()
//                .map(r -> r.getFeature(SecStrucState.class, FeatureNames.SECONDARY_STRUCTURE_STATES))
//                .forEach(System.out::println);
    }

    private void detectStrands() {
        // Find all the beta bridges of the structure
        findBridges();
        // Create Ladders
        createLadders();
        // Detect beta bulges between ladders
        connectLadders();
        // AND store SS assignments for Sheets, Strands and Bridges
        updateSheets();
    }

    private SecStrucState getState(Residue residue) {
        return residue.getFeature(SecStrucState.class, FeatureNames.SECONDARY_STRUCTURE_STATES);
    }

    private void updateSheets() {
        for(Ladder ladder : ladders){
            for (int lcount = ladder.getFrom(); lcount <= ladder.getTo(); lcount++) {
                SecStrucState state = getState(residues.get(lcount));
                DSSPSecondaryStructureElement stype = state.getSecondaryStructure();

                int diff = ladder.getFrom() - lcount;
                int l2count = ladder.getLfrom() - diff ;

                SecStrucState state2 = getState(residues.get(l2count));
                DSSPSecondaryStructureElement stype2 = state2.getSecondaryStructure();

                if(ladder.getFrom() != ladder.getTo()) {
                    setSecStrucType(lcount, DSSPSecondaryStructureElement.EXTENDED);
                    setSecStrucType(l2count, DSSPSecondaryStructureElement.EXTENDED);
                } else {
                    if(!stype.isHelixType() && (!stype.equals(DSSPSecondaryStructureElement.EXTENDED)))
                        setSecStrucType(lcount, DSSPSecondaryStructureElement.BRIDGE);

                    if(!stype2.isHelixType() && (!stype2.equals(DSSPSecondaryStructureElement.EXTENDED)))
                        setSecStrucType(l2count, DSSPSecondaryStructureElement.BRIDGE);
                }
            }

            // Check if two ladders are connected. both sides are 'E'

            if(ladder.getConnectedTo() == 0) {
                continue;
            }
            Ladder conladder = ladders.get(ladder.getConnectedTo());

            if (ladder.getBtype().equals(BridgeType.ANTIPARALLEL)) {
				/* set one side */
                for(int lcount = ladder.getFrom(); lcount <= conladder.getTo(); lcount++) {
                    setSecStrucType(lcount, DSSPSecondaryStructureElement.EXTENDED);
                }
				/* set other side */
                for (int lcount = conladder.getLto(); lcount <= ladder.getLfrom(); lcount++) {
                    setSecStrucType(lcount, DSSPSecondaryStructureElement.EXTENDED);
                }

            } else {
				/* set one side */
                for(int lcount = ladder.getFrom(); lcount <= conladder.getTo(); lcount++) {
                    setSecStrucType(lcount, DSSPSecondaryStructureElement.EXTENDED);
                }
				/* set other side */
                for(int lcount = ladder.getLfrom(); lcount <= conladder.getLto(); lcount++) {
                    setSecStrucType(lcount, DSSPSecondaryStructureElement.EXTENDED);
                }
            }
        }
    }

    private void connectLadders() {
        for(int i = 0; i < ladders.size(); i++) {
            for(int j = i; j < ladders.size(); j++) {
                Ladder l1 = ladders.get(i);
                Ladder l2 = ladders.get(j);
                System.out.println(l1 + " : " + l2);
                if (hasBulge(l1, l2)) {
                    l1.setConnectedTo(i);
                    l2.setConnectedFrom(j);
                    logger.debug("Bulge for: ({}, {})", i, j);
                }
            }
        }
    }

    /**
     * For beta structures, we define explicitly: a bulge-linked
     * ladder consists of two (perfect) ladder or bridges of the
     * same type connected by at most one extra getResidue on one
     * strand and at most four extra residues on the other strand,
     * all residues in bulge-linked ladders are marked "E,"
     * including the extra residues.
     */
    private boolean hasBulge(Ladder l1, Ladder l2) {
        boolean bulge = ((l1.getBtype().equals(l2.getBtype())) &&
                (l2.getFrom() - l1.getTo() < 6) &&
                (l1.getTo() < l2.getFrom()) &&
                (l2.getConnectedTo() == 0));

        if(!bulge) {
            return bulge;
        }

        switch(l1.getBtype()){
            case PARALLEL:
                return ((l2.getLfrom() - l1.getLto() > 0) &&
                        (((l2.getLfrom() - l1.getLto() < 6) &&
                                (l2.getFrom() - l1.getTo() < 3)) ||
                                (l2.getLfrom() - l1.getLto() < 3)));
            case ANTIPARALLEL:
                return ((l1.getLfrom() - l2.getLto() > 0) &&
                        (((l1.getLfrom() - l2.getLto() < 6) &&
                                (l2.getFrom() - l1.getTo() < 3)) ||
                                (l1.getLfrom() - l2.getLto() < 3)));
            default: throw new IllegalArgumentException("case " + l1.getBtype() + " not supported");
        }
    }

    private void createLadders(){
        for (BetaBridge b : bridges) {
            boolean found = false;
            for(Ladder ladder : ladders) {
                if(shouldExtendLadder(ladder, b)) {
                    found = true;
                    ladder.setTo(ladder.getTo() + 1); //we go forward in this direction
                    switch(b.getType()){
                        case PARALLEL:
                            ladder.setLto(ladder.getLto() + 1); //increment second strand
                            break;
                        case ANTIPARALLEL:
                            ladder.setLfrom(ladder.getLfrom() - 1); //decrement second strand
                            break;
                    }
                    break;
                }
            }
            if (!found) {
                //Create new ladder with a single Bridge
                Ladder l = new Ladder();
                l.setFrom(b.getPartner1());
                l.setTo(b.getPartner1());
                l.setLfrom(b.getPartner2());
                l.setLto(b.getPartner2());
                l.setBtype(b.getType());
                ladders.add(l);
            }
        }
    }

    /**
     * Conditions to extend a ladder with a given beta Bridge:
     * <li>The bridge and ladder are of the same type.
     * <li>The smallest bridge getResidue is sequential to the first
     * 		strand ladder.
     * <li>The second bridge getResidue is either sequential (parallel)
     * 		or previous (antiparallel) to the second strand of the ladder
     * </li>
     * @param ladder the ladder candidate to extend
     * @param b the beta bridge that would extend the ladder
     * @return true if the bridge b extends the ladder
     */
    private boolean shouldExtendLadder(Ladder ladder, BetaBridge b) {
        //Only extend if they are of the same type
        boolean sameType = b.getType().equals(ladder.getBtype());
        if (!sameType) {
            return false;
        }

        //Only extend if getResidue 1 is sequential to ladder strand
        boolean sequential = (b.getPartner1() == ladder.getTo() + 1);
        if (!sequential) {
            return false;
        }

        switch(b.getType()){
            case PARALLEL:
                //Residue 2 should be sequential to second strand
                if (b.getPartner2() == ladder.getLto() + 1) {
                    return true;
                }
            case ANTIPARALLEL:
                //Residue 2 should be previous to second strand
                if (b.getPartner2() == ladder.getLfrom() - 1) {
                    return true;
                }
            default: return false;
        }
    }

    /**
     * Two nonoverlapping stretches of three residues each, i-1,i,i+1 and
     * j-1,j,j+1, form either a parallel or antiparallel bridge, depending on
     * which of two basic patterns is matched. We assign a bridge between
     * residues i and j if there are two H bonds characteristic of beta-
     * structure; in particular:
     * <p>
     * Parallel Bridge(i,j) =: [Hbond(i-1,j) and Hbond(j,i+1)]
     * 							or [Hbond(j-1,i) and Hbond(i,j+1)]
     * <p>
     * Antiparallel Bridge(i,j) =: [Hbond(i,j) and Hbond(j,i)]
     * 								or [Hbond(i-1,j+1) and Hbond(j-1,i+1)]
     *
     * Optimised to use the contact set
     */
    private void findBridges() {
        List<int[]> outList = new ArrayList<>();

        for (int i = 0; i < residues.size() - 1; i++) {
            for (int j = i + 1; j < residues.size(); j++) {
                Residue res1 = residues.get(i);
                Residue res2 = residues.get(j);
                double distance = LinearAlgebra3D.distanceFast(
                        res1.getAlphaCarbon().getCoordinates(),
                        res2.getAlphaCarbon().getCoordinates());
                if (distance > CA_MIN_DIST) {
                    continue;
                }
                // If i>j switch them over
                if(i > j){
                    // Switch them over
                    int old = i;
                    i = j;
                    j = old;
                }
                // Only these
                if(j < i + 3){
                    continue;
                }
                // If it's the first
                if(i == 0){
                    continue;
                }
                if(j == 0){
                    continue;
                }
                // If it's the last
                if(i == residues.size() - 1){
                    continue;
                }
                if(j == residues.size() - 1){
                    continue;
                }

                int[] thisPair = new int[] { i, j };
                outList.add(thisPair);
            }
        }

        Collections.sort(outList, (o1, o2) -> {
            if(o1[0] < o2[0]) {
                return -1;
            } else if(o1[0] > o2[0]) {
                return +1;
            } else {
                if(o1[1] < o2[1]) {
                    return -1;
                } else if(o1[1] > o2[1]) {
                    return +1;
                } else {
                    return 0;
                }
            }
        });

        for(int[] p : outList){
            int i = p[0];
            int j = p[1];
            BridgeType btype = null;
            // Now do the bonding
            if((isBonded(i-1,j) && isBonded(j,i+1)) ||
                    (isBonded(j-1,i) && isBonded(i,j+1))) {
                btype = BridgeType.PARALLEL;
            }
            else if ((isBonded(i,j) && isBonded(j,i)) ||
                    (isBonded(i-1,j+1) && (isBonded(j-1,i+1)))) {
                btype = BridgeType.ANTIPARALLEL;
            }
            if (btype != null){
                registerBridge(i, j, btype);
            }
        }
    }

    private void registerBridge(int i, int j, BridgeType btype) {
        BetaBridge bridge = new BetaBridge(i, j, btype);

        boolean b1 = getState(residues.get(i)).addBridge(bridge);
        boolean b2 = getState(residues.get(j)).addBridge(bridge);

        if (!b1 && !b2) {
            logger.warn("Ignoring Bridge between residues {} and {}. DSSP assignment might differ.", i, j);
        }

        bridges.add(bridge);
    }


    private void detectBends() {
        f1: for (int i = 2; i < residues.size() - 2; i++) {
            // Check if all atoms form peptide bonds (backbone discontinuity)
            for (int k = 0; k < 4; k++) {
                int index = i + k - 2;
                Atom c = residues.get(index).getBackboneCarbon();
                Atom n = residues.get(index + 1).getBackboneNitrogen();
                // Peptide bond C-N
                if (LinearAlgebra3D.distanceFast(c.getCoordinates(), n.getCoordinates()) >
                        MAX_PEPTIDE_BOND_LENGTH_SQUARED) {
                    continue f1;
                }
            }

            Atom caim2 = residues.get(i - 2).getAlphaCarbon();
            Atom cag = residues.get(i).getAlphaCarbon();
            Atom caip2 = residues.get(i + 2).getAlphaCarbon();

            // Create vectors ( Ca i to Ca i-2 ) ; ( Ca i to CA i + 2 )
            double[] caminus2 = LinearAlgebra3D.subtract(caim2.getCoordinates(), cag.getCoordinates());
            double[] caplus2 = LinearAlgebra3D.subtract(cag.getCoordinates(), caip2.getCoordinates());

            double angle = LinearAlgebra3D.angle(caminus2, caplus2);

            SecStrucState state = getState(residues.get(i));
            state.setKappa((float) angle);

            // Angles = 360 should be discarded
            if (angle > 70.0 && angle < 359.99) {
                setSecStrucType(i, DSSPSecondaryStructureElement.BEND);
                state.setBend(true);
            }
        }
    }

    private void buildHelices() {
        // Alpha-helix (i+4), 3-10-helix (i+3), Pi-helix (i+5)
        checkSetHelix(4, DSSPSecondaryStructureElement.ALPHA_HELIX);
        checkSetHelix(3, DSSPSecondaryStructureElement.THREE10HELIX);
        checkSetHelix(5, DSSPSecondaryStructureElement.PIHELIX);

        checkSetTurns();
    }

    private void checkSetTurns() {
        DSSPSecondaryStructureElement type = DSSPSecondaryStructureElement.TURN;

        for(int idx = 0; idx < 3; idx++) {
            for(int i = 0; i < residues.size() - 1; i++) {
                SecStrucState state = getState(residues.get(i));
                char[] turn = state.getTurn();

                // Any turn opening matters
                if(turn[idx] == '>' || turn[idx] == 'X') {
                    // Mark following n residues as turn
                    for(int k = 1; k < idx + 3; k++) {
                        setSecStrucType(i + k, type);
                    }
                    if(!DSSP_HELICES) {
                        setSecStrucType(i, type);
                        setSecStrucType(i + idx + 3, type);
                    }
                }
            }
        }
    }

    /**
     * Set the new type only if it has more preference than the current getResidue
     * SS type.
     *
     * @param pos the getResidue index to manipulate
     * @param type the type to assign
     */
    private void setSecStrucType(int pos, DSSPSecondaryStructureElement type) {
        SecStrucState ss = getState(residues.get(pos));
        // more favorable according to DSSP ranking
        if (type.compareTo(ss.getSecondaryStructure()) > 0) {
            ss.setSecondaryStructure(type);
        }
    }

    /**
     * A minimal helix is defined by two consecutive n-turns. For example, a
     * 4-helix, of minimal length 4 from residues i to (i+3), requires turns (of
     * type 4) at residues (i-1) and i.
     * <p>
     * Note that the orignal DSSP implementation does not assign helix type to
     * getResidue (i-1) and getResidue (i+n+1), although they contain a helix turn. As
     * they state in the original paper, "the helices are one getResidue shorter
     * than they would be according to rule 6.3 of IUPAC-IUB".
     *
     * @param n
     * @param type
     */
    private void checkSetHelix(int n, DSSPSecondaryStructureElement type) {
        int idx = n - 3;
        logger.debug("Set helix {} {} {}", type, n, idx);

        for (int i = 1; i < residues.size() - n; i++) {
            SecStrucState state = getState(residues.get(i));
            SecStrucState previousState = getState(residues.get(i - 1));

            // Check that no other helix was assgined to this range
            if (state.getSecondaryStructure().compareTo(type) > 0) {
                continue;
            }
            if (getState(residues.get(i + 1)).getSecondaryStructure().compareTo(type) > 0) {
                continue;
            }

            char turn = state.getTurn()[idx];
            char pturn = previousState.getTurn()[idx];

            // Two consecutive n-turns present to define a n-helix
            if ((turn == '>' || turn == 'X') && (pturn == '>' || pturn == 'X')) {
                // Mark following n residues as turn
                for (int k = 0; k < n; k++) {
                    setSecStrucType(i + k, type);
                }
                if (!DSSP_HELICES) {
                    setSecStrucType(i - 1, type);
                    setSecStrucType(i + n, type);
                }
            }
        }
    }

    /**
     * Detect helical turn patterns.
     */
    private void calculateTurns() {
        for (int i = 0; i < residues.size(); i++) {
            Residue residue1 = residues.get(i);
            for (int turn = 3; turn <= 5; turn++) {
                if (i + turn >= residues.size()) {
                    continue;
                }

                Residue residue2 = residues.get(i + turn);
                // Check for H bond from NH(i+n) to CO(i)
                if (isBonded(i, i + turn)) {
                    logger.debug("Turn at ({}, {}, {}", i, i + turn, turn);
                    getState(residue1).setTurn('>', turn);
                    getState(residue2).setTurn('<', turn);
                    // Bracketed residues get the helix number
                    for (int j = i + 1; j < i + turn; j++) {
                        int t = turn;
                        char helix = String.valueOf(t).charAt(0);
                        getState(residues.get(j)).setTurn(helix, turn);
                    }
                }
            }
        }
    }

    /**
     * Test if two groups are forming an H-Bond. The bond tested is from the CO
     * of group i to the NH of group j. Acceptor (i) and donor (j). The donor of
     * i has to be j, and the acceptor of j has to be i. DSSP defines H-Bonds if
     * the energy < -500 cal/mol.
     *
     * @param i
     *            group one
     * @param j
     *            group two
     * @return flag if the two are forming an Hbond
     */
    private boolean isBonded(int i, int j) {
        Residue residue1 = residues.get(i);
        Residue residue2 = residues.get(j);
        SecStrucState state1 = getState(residue1);
        SecStrucState state2 = getState(residue2);

        double don1e = state1.getDonor1().getEnergy();
        double don2e = state1.getDonor2().getEnergy();
        double acc1e = state2.getAccept1().getEnergy();
        double acc2e = state2.getAccept2().getEnergy();

        Residue don1p = state1.getDonor1().getPartner();
        Residue don2p = state1.getDonor2().getPartner();
        Residue acc1p = state2.getAccept1().getPartner();
        Residue acc2p = state2.getAccept2().getPartner();

        // Either donor from i is j, or accept from j is i
        boolean hbond = (residue2.equals(don1p) && don1e < HBONDHIGHENERGY) || (residue2.equals(don2p) && don2e < HBONDHIGHENERGY)
                || (residue1.equals(acc1p) && acc1e < HBONDHIGHENERGY) || (residue1.equals(acc2p) && acc2e < HBONDHIGHENERGY);

        if (hbond) {
            logger.debug("*** H-bond from CO of {} to NH of {}", i, j);
            return true;
        }
        return false;
    }

    private void calculateDihedralAngles() {
        // dihedral angles
        // phi: C-N-CA-C
        // psi: N-CA-C-N
        // Chi1: N-CA-CB-CG, N-CA-CB-OG(SER),N-CA-CB-OG1(Thr),
        // N-CA-CB-CG1(ILE/VAL), N-CA-CB-SG(CYS)
        // Omega: CA-C-N-CA
        for (int i = 0; i < residues.size() - 1; i++) {
            Residue res1 = residues.get(i);
            Residue res2 = residues.get(i + 1);

            Atom n1 = res1.getBackboneNitrogen();
            Atom ca1 = res1.getAlphaCarbon();
            Atom c1 = res1.getBackboneCarbon();

            Atom n2 = res2.getBackboneNitrogen();
            Atom ca2 = res2.getAlphaCarbon();
            Atom c2 = res2.getBackboneCarbon();

            double phi = LinearAlgebra3D.torsionAngle(c1.getCoordinates(), n2.getCoordinates(), ca2.getCoordinates(), c2.getCoordinates());
            double psi = LinearAlgebra3D.torsionAngle(n1.getCoordinates(), ca1.getCoordinates(), c1.getCoordinates(), n2.getCoordinates());
            double omega = LinearAlgebra3D.torsionAngle(ca1.getCoordinates(), c1.getCoordinates(), n2.getCoordinates(), ca2.getCoordinates());

            SecStrucState state1 = getState(res1);
            SecStrucState state2 = getState(res2);

            state2.setPhi(phi);
            state1.setPsi(psi);
            state1.setOmega(omega);
        }
    }

    /**
     * Calculate the HBonds between different groups. see Creighton page 147 f
     * Modified to use only the contact map
     */
    private void calculateHBonds() {
//        residues.
        double squaredDistanceCutoff = CA_MIN_DIST * CA_MIN_DIST;
        for (int i = 0; i < residues.size() - 1; i++) {
            for (int j = i + 1; j < residues.size(); j++) {
                Residue res1 = residues.get(i);
                Residue res2 = residues.get(j);
                double squaredDistance = LinearAlgebra3D.distanceFast(
                        res1.getAlphaCarbon().getCoordinates(),
                        res2.getAlphaCarbon().getCoordinates());
                if (squaredDistance > squaredDistanceCutoff) {
                    continue;
                }
                checkAddHBond(res1, res2);
                if (j != (i + 1)) {
                    checkAddHBond(res2, res1);
                }
            }
        }
    }

    private void checkAddHBond(Residue residue1, Residue residue2) {
        if (residue1.getAminoAcid().equals(AminoAcid.PROLINE)) {
            logger.debug("Ignore: PRO {}", residue1.getResidueNumber());
            return;
        }
        if (lacksBackboneHydrogen(residue1)) {
            logger.debug("Residue {} has no H", residue1.getResidueNumber());
            return;
        }

        Pair<Residue, Residue> residuePair = new Pair<>(residue1, residue2);
        trackHBondEnergy(residuePair, calculateHBondEnergy(residuePair));
    }

    /**
     * Store Hbonds in the Groups. DSSP allows two HBonds per amino acid to
     * allow bifurcated bonds.
     */
    private void trackHBondEnergy(Pair<Residue, Residue> residuePair, double energy) {
        Residue res1 = residuePair.getFirst();
        Residue res2 = residuePair.getSecond();
        if (res1.getAminoAcid().equals(AminoAcid.PROLINE)) {
            logger.debug("Ignore: PRO {}", res1.getResidueNumber());
            return;
        }

        // try to get entries or create new ones with coil secondary structure
        SecStrucState state1 = getState(res1);
        SecStrucState state2 = getState(res2);

        double acc1e = state1.getAccept1().getEnergy();
        double acc2e = state1.getAccept2().getEnergy();

        double don1e = state2.getDonor1().getEnergy();
        double don2e = state2.getDonor2().getEnergy();

        // Acceptor: N-H-->O
        if (energy < acc1e) {
            logger.debug("{} < {}", energy, acc1e);
            state1.setAccept2(state1.getAccept1());

            HBond bond = new HBond();
            bond.setEnergy(energy);
            bond.setPartner(res2);

            state1.setAccept1(bond);
        } else if (energy < acc2e) {
            logger.debug("{} < {}", energy, acc2e);

            HBond bond = new HBond();
            bond.setEnergy(energy);
            bond.setPartner(res2);

            state1.setAccept2(bond);
        }

        // The other side of the bond: donor O-->N-H
        if (energy < don1e) {
            logger.debug("{} < {}", energy, don1e);
            state2.setDonor2(state2.getDonor1());

            HBond bond = new HBond();
            bond.setEnergy(energy);
            bond.setPartner(res1);

            state2.setDonor1(bond);
        } else if (energy < don2e) {
            logger.debug("{} < {}", energy, don2e);

            HBond bond = new HBond();
            bond.setEnergy(energy);
            bond.setPartner(res1);

            state2.setDonor2(bond);
        }
    }

    /**
     * Calculate HBond energy of two groups in cal/mol see Creighton page 147 f
     * <p>
     * Jeffrey, George A., An introduction to hydrogen bonding, Oxford
     * University Press, 1997. categorizes hbonds with donor-acceptor distances
     * of 2.2-2.5 &aring; as "strong, mostly covalent", 2.5-3.2 &aring; as
     * "moderate, mostly electrostatic", 3.2-4.0 &aring; as "weak,
     * electrostatic". Energies are given as 40-14, 15-4, and <4 kcal/mol
     * respectively.
     * @param residuePair
     * @return the energy of this h-bond
     */
    private double calculateHBondEnergy(Pair<Residue, Residue> residuePair) {
        Residue res1 = residuePair.getFirst();
        Residue res2 = residuePair.getSecond();
        Atom nAtom = res1.getBackboneNitrogen();
        double[] n = nAtom.getCoordinates();
        double[] h = res1.getBackboneHydrogen().getCoordinates();

        Atom oAtom = res2.getBackboneOxygen();
        double[] o = oAtom.getCoordinates();
        double[] c = res2.getBackboneCarbon().getCoordinates();

        double dno = LinearAlgebra3D.distance(o, n);
        double dhc = LinearAlgebra3D.distance(c, h);
        double dho = LinearAlgebra3D.distance(o, h);
        double dnc = LinearAlgebra3D.distance(c, n);

//        logger.debug("     cccc: " + res1.getResidueNumber() + " " + res1.getAminoAcid() + " " +
//                res2.getResidueNumber() + " " + res2.getAminoAcid() + String.format( " O (" + oAtom.getPdbSerial() +
//                ")..N (" + nAtom.getPdbSerial() + "):%4.1f  |  ho:%4.1f - hc:%4.1f + nc:%4.1f - no:%4.1f ", dno, dho,
//                dhc, dnc, dno));

        // there seems to be a contact!
        if ((dno < MINDIST) || (dhc < MINDIST) || (dnc < MINDIST) || (dno < MINDIST)) {
            return HBONDLOWENERGY;
        }

        double e1 = Q / dho - Q / dhc;
        double e2 = Q / dnc - Q / dno;

        double energy = e1 + e2;

//        logger.debug(String.format("      N (%d) O(%d): %4.1f : %4.2f ", nAtom.getPdbSerial(),
//                oAtom.getPdbSerial(), (float) dno, energy));

        // Avoid too strong energy
        if (energy > HBONDLOWENERGY) {
            return energy;
        }

        return HBONDLOWENERGY;
    }

    /**
     * Calculate the coordinates of the H atoms. They are usually missing in the
     * PDB files as only few experimental methods allow to resolve their
     * location.
     */
    private void calculateHAtoms() {
        Fragment.fragmentsOf(residues, 2).filter(this::lacksBackboneHydrogen).forEach(this::calcSimpleH);
    }

    /**
     * <code>true</code> iff the 2nd getResidue of this fragment lacks the backbone hydrogen
     * @param fragmentOfSize2 2 consecutively connected residues
     * @return true if no backbone hydrogen is present
     */
    private boolean lacksBackboneHydrogen(Fragment<Residue> fragmentOfSize2) {
        return lacksBackboneHydrogen(fragmentOfSize2.getElement(1));
    }

    /**
     * <code>true</code> iff the 2nd getResidue of this fragment lacks the backbone hydrogen
     * @param residue the getResidue to check
     * @return true if no backbone hydrogen is present
     */
    private boolean lacksBackboneHydrogen(Residue residue) {
        return !residue.findBackboneHydrogen().isPresent();
    }

    /**
     * Places a pseudo hydrogen atoms between 2 residues, when it is lacking.
     * @param fragmentOfSize2 the fragment where the atom shall be placed
     */
    private void calcSimpleH(Fragment<Residue> fragmentOfSize2) {
        try {
            double[] c = fragmentOfSize2.getElement(0).getBackboneCarbon().getCoordinates();
            double[] o = fragmentOfSize2.getElement(0).getBackboneOxygen().getCoordinates();
            double[] n = fragmentOfSize2.getElement(1).getBackboneNitrogen().getCoordinates();

            double[] xyz = LinearAlgebra3D.add(n, LinearAlgebra3D.divide(LinearAlgebra3D.subtract(c, o),
                    LinearAlgebra3D.distance(o, c)));

            // pdbSerial of minimal int value flags them as pseudo-hydrogen
            fragmentOfSize2.getElement(1).addAtom(new Atom("H", Element.H, xyz));
        } catch (NoSuchElementException e) {
            logger.warn(e.getLocalizedMessage());
        }
    }
}