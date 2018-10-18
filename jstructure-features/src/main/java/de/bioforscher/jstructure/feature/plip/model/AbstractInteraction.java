package de.bioforscher.jstructure.feature.plip.model;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.selection.SelectionException;

import java.util.List;

public abstract class AbstractInteraction implements Interaction {
    private final List<Atom> interactingAtoms1;
    private final List<Atom> interactingAtoms2;
    private final Group interactingGroup1;
    private final Group interactingGroup2;
    private final boolean polymerInteraction;
    private final boolean ligandInteraction;

    AbstractInteraction(List<Atom> interactingAtoms1,
                        List<Atom> interactingAtoms2,
                        Group interactingGroup1,
                        Group interactingGroup2) {
        this.interactingAtoms1 = interactingAtoms1;
        this.interactingGroup1 = interactingGroup1;
        this.interactingAtoms2 = interactingAtoms2;
        this.interactingGroup2 = interactingGroup2;
        this.polymerInteraction = !interactingGroup1.isLigand() && !interactingGroup2.isLigand();
        this.ligandInteraction = !polymerInteraction;
    }

    @Override
    public boolean isPolymerInteraction() {
        return polymerInteraction;
    }

    @Override
    public boolean isLigandInteraction() {
        return ligandInteraction;
    }

    @Override
    public boolean isDirectedInteraction() {
        return false;
    }

    @Override
    public List<Atom> getInteractingAtoms1() {
        return interactingAtoms1;
    }

    @Override
    public List<Atom> getInteractingAtoms2() {
        return interactingAtoms2;
    }

    public Group getInteractingGroup1() {
        return interactingGroup1;
    }

    public Group getInteractingGroup2() {
        return interactingGroup2;
    }

    @Override
    public boolean containsGroup(Group group) {
        return interactingGroup1.equals(group) || interactingGroup2.equals(group);
    }

    @Override
    public Group getOpposingGroup(Group group) {
        if(interactingGroup1.equals(group)) {
            return interactingGroup2;
        }
        if(interactingGroup2.equals(group)) {
            return interactingGroup1;
        }
        throw new SelectionException("interaction " + this + " does not contain group " + group);
    }

    @Override
    public String toString() {
        return getClass().getSimpleName() + " " + interactingGroup1 + " <-> " + interactingGroup2 + " : " + interactingAtoms1 + " <-> " + interactingAtoms2;
    }
}
