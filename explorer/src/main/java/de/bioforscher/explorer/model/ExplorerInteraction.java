package de.bioforscher.explorer.model;

import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.parser.plip.interaction.*;

/**
 * The reduced representation of a PLIP interaction.
 * Created by bittrich on 2/22/17.
 */
@SuppressWarnings("unused")
public class ExplorerInteraction {
    private String a1, a2;

    public ExplorerInteraction() {

    }

    ExplorerInteraction(HalogenBond halogenBond) {
        this.a1 = convertAtomToPVSelection(halogenBond.getAcceptor());
        this.a2 = convertAtomToPVSelection(halogenBond.getDonor());
    }

    ExplorerInteraction(HydrogenBond hydrogenBond) {
        this.a1 = convertAtomToPVSelection(hydrogenBond.getAcceptor());
        this.a2 = convertAtomToPVSelection(hydrogenBond.getDonor());
    }

    ExplorerInteraction(HydrophobicInteraction hydrophobicInteraction) {
        this.a1 = convertAtomToPVSelection(hydrophobicInteraction.getAtom1());
        this.a2 = convertAtomToPVSelection(hydrophobicInteraction.getAtom2());
    }

    ExplorerInteraction(MetalComplex metalComplex) {
        this.a1 = convertAtomToPVSelection(metalComplex.getAtom1());
        this.a2 = convertAtomToPVSelection(metalComplex.getAtom2());
    }

    ExplorerInteraction(Pair<Atom, Atom> atoms) {
        this.a1 = convertAtomToPVSelection(atoms.getLeft());
        this.a2 = convertAtomToPVSelection(atoms.getRight());
    }

    ExplorerInteraction(WaterBridge waterBridge) {
        this.a1 = convertAtomToPVSelection(waterBridge.getAtom1());
        this.a2 = convertAtomToPVSelection(waterBridge.getAtom2());
    }

    private String convertAtomToPVSelection(Atom atom) {
        return atom.getParentGroup().getParentChain().getChainId() + "." + atom.getParentGroup().getResidueNumber() + "." + atom.getName();
    }

    public String getA1() {
        return a1;
    }

    public String getA2() {
        return a2;
    }
}
