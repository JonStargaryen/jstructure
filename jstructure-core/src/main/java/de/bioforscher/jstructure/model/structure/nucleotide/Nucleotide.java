package de.bioforscher.jstructure.model.structure.nucleotide;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.GroupPrototype;
import de.bioforscher.jstructure.model.structure.ResidueNumber;

/**
 * Created by bittrich on 5/30/17.
 */
public abstract class Nucleotide extends Group implements StandardNucleotideIndicator {
    private Atom op3;
    private Atom p;
    private Atom op1;
    private Atom op2;
    private Atom o5prime;
    private Atom c5prime;
    private Atom c4prime;
    private Atom o4prime;
    private Atom c3prime;
    private Atom o3prime;
    private Atom c2prime;
    private Atom c1prime;
    private Atom n1;
    private Atom c2;
    private Atom n3;
    private Atom c4;
    private Atom c5;
    private Atom c6;

    //TODO grouping of atom names & family enum

    Nucleotide(Nucleotide nucleotide) {
        super(nucleotide);
        atoms().forEach(this::addAtomInternal);
    }

    Nucleotide(String threeLetterCode,
              ResidueNumber residueNumber,
              boolean ligand) {
        super(threeLetterCode,
                residueNumber,
                ligand);
    }

    Nucleotide(String threeLetterCode,
              ResidueNumber residueNumber) {
        this(threeLetterCode, residueNumber, false);
    }

    Nucleotide(GroupPrototype groupPrototype,
              ResidueNumber residueNumber,
              boolean ligand) {
        super(groupPrototype,
                residueNumber,
                ligand);
    }

    Nucleotide(GroupPrototype groupPrototype,
              ResidueNumber residueNumber) {
        this(groupPrototype, residueNumber, false);
    }

    @Override
    protected void addAtomInternal(Atom atom) {
        if(atom.getName().equals("OP3") && op3 == null) {
            op3 = atom;
        }
        if(atom.getName().equals("P") && p == null) {
            p = atom;
        }
        if(atom.getName().equals("OP1") && op1 == null) {
            op1 = atom;
        }
        if(atom.getName().equals("OP2") && op2 == null) {
            op2 = atom;
        }
        if(atom.getName().equals("\"O5'\"") && o5prime == null) {
            o5prime = atom;
        }
        if(atom.getName().equals("\"C5'\"") && c5prime == null) {
            c5prime = atom;
        }
        if(atom.getName().equals("\"C4'\"") && c4prime == null) {
            c4prime = atom;
        }
        if(atom.getName().equals("\"O4'\"") && o4prime == null) {
            o4prime = atom;
        }
        if(atom.getName().equals("\"C3'\"") && c3prime == null) {
            c3prime = atom;
        }
        if(atom.getName().equals("\"O3'\"") && o3prime == null) {
            o3prime = atom;
        }
        if(atom.getName().equals("\"C2'\"") && c2prime == null) {
            c2prime = atom;
        }
        if(atom.getName().equals("\"C1'\"") && c1prime == null) {
            c1prime = atom;
        }
        if(atom.getName().equals("N1") && n1 == null) {
            n1 = atom;
        }
        if(atom.getName().equals("C2") && c2 == null) {
            c2 = atom;
        }
        if(atom.getName().equals("N3") && n3 == null) {
            n3 = atom;
        }
        if(atom.getName().equals("C4") && c4 == null) {
            c4 = atom;
        }
        if(atom.getName().equals("C5") && c5 == null) {
            c5 = atom;
        }
        if(atom.getName().equals("C6") && c6 == null) {
            c6 = atom;
        }
        addBaseAtom(atom);
    }

    protected abstract void addBaseAtom(Atom atom);
}
