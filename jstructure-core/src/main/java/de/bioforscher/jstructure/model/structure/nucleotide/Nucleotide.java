package de.bioforscher.jstructure.model.structure.nucleotide;

import de.bioforscher.jstructure.model.structure.*;
import de.bioforscher.jstructure.model.structure.identifier.ResidueIdentifier;

import java.util.Optional;
import java.util.stream.Stream;

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
              ResidueIdentifier residueIdentifier,
              boolean ligand) {
        super(threeLetterCode,
                residueIdentifier,
                ligand);
    }

    Nucleotide(String threeLetterCode,
              ResidueIdentifier residueIdentifier) {
        this(threeLetterCode, residueIdentifier, false);
    }

    Nucleotide(GroupPrototype groupPrototype,
              ResidueIdentifier residueIdentifier,
              boolean ligand) {
        super(groupPrototype,
                residueIdentifier,
                ligand);
    }

    Nucleotide(GroupPrototype groupPrototype,
              ResidueIdentifier residueIdentifier) {
        this(groupPrototype, residueIdentifier, false);
    }

    public enum Family implements GroupFamily {
        ADENOSINE(Adenosine.class,
                Adenosine.GROUP_PROTOTYPE),
        CYTIDINE(Cytidine.class,
                Cytidine.GROUP_PROTOTYPE),
        DEOXYADENOSINE(Deoxyadenosine.class,
                Deoxyadenosine.GROUP_PROTOTYPE),
        DEOXYCYTIDINE(Deoxycytidine.class,
                Deoxycytidine.GROUP_PROTOTYPE),
        DEOXYGUANOSINE(Deoxyguanosine.class,
                Deoxyguanosine.GROUP_PROTOTYPE),
        GUANOSINE(Guanosine.class,
                Guanosine.GROUP_PROTOTYPE),
        THYMIDINE(Thymidine.class,
                Thymidine.GROUP_PROTOTYPE),
        URIDINE(Uridine.class,
                Uridine.GROUP_PROTOTYPE),
        UNKNOWN_NUCLEOTIDE(UnknownNucleotide.class,
                UnknownNucleotide.GROUP_PROTOTYPE);

        private Class<? extends Nucleotide> representingClass;
        private GroupPrototype groupPrototype;

        Family(Class<? extends Nucleotide> representingClass,
               GroupPrototype groupPrototype) {
            this.representingClass = representingClass;
            this.groupPrototype = groupPrototype;
        }

        @Override
        public Class<? extends Nucleotide> getRepresentingClass() {
            return representingClass;
        }

        @Override
        public GroupPrototype getGroupPrototype() {
            return groupPrototype;
        }

        public static Nucleotide createNucleotide(String pdbName, ResidueIdentifier residueIdentifier, boolean ligand) {
            Class<? extends Nucleotide> representingClass = resolveThreeLetterCode(pdbName).representingClass;
            // use special constructor for UnknownNucleotide
            if(representingClass.isAssignableFrom(UnknownNucleotide.class)) {
                return new UnknownNucleotide(pdbName, residueIdentifier, ligand);
            } else {
                try {
                    return representingClass.getConstructor(ResidueIdentifier.class, boolean.class).newInstance(residueIdentifier, ligand);
                } catch (Exception e) {
                    throw new RuntimeException("creation of Nucleotide instance failed", e);
                }
            }
        }

        public static Nucleotide.Family resolveThreeLetterCode(String threeLetterCode) {
            return Stream.of(Nucleotide.Family.values())
                    .filter(nucleotide -> threeLetterCode.equalsIgnoreCase(nucleotide.getThreeLetterCode()))
                    .findFirst()
                    .orElse(Nucleotide.Family.UNKNOWN_NUCLEOTIDE);
        }

        public static Nucleotide.Family resolveGroupPrototype(GroupPrototype groupPrototype) {
            GroupPrototype.PolymerType polymerType = groupPrototype.getPolymerType();
            if(groupPrototype.getPolymerType() != GroupPrototype.PolymerType.NA_LINKING) {
                throw new UnsupportedOperationException("method only supported for nucleotides - group '" +
                        groupPrototype.getThreeLetterCode() + "' has polymer type: " + polymerType);
            }

            return Stream.of(Nucleotide.Family.values())
                    .filter(nucleotide -> groupPrototype.equals(nucleotide.getGroupPrototype()))
                    .findFirst()
                    .orElse(UNKNOWN_NUCLEOTIDE);
        }

        //TODO retrieval functions for DNA-/RNA-bases etc

        @Override
        public String getThreeLetterCode() {
            return groupPrototype.getThreeLetterCode();
        }
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

    /**
     * Access to the previous nucleotide in the chain.
     * <code>getPreviousNucleotide() is equivalent to
     * </code><blockquote><pre>
     * getNucleotideWithOffset(-1)
     * </pre></blockquote>
     * @see #getNucleotideWithOffset(int)
     * @return an optional wrapping the previous nucleotide
     */
    public Optional<Nucleotide> getPreviousNucleotide() {
        return getNucleotideWithOffset(-1);
    }

    /**
     * Access to the nucleotide with a custom offset in the chain. Zero will return the current instance, negative
     * values will navigate towards the N-terminus, positive values to the C-terminus. Wrapped as optional as the
     * operation is by no means guaranteed to succeed, the group could not exist (as the current instance is a terminal
     * nucleotide) or be no {@link Nucleotide}.
     * @see #getPreviousNucleotide()
     * @see #getNextNucleotide()
     * @return an optional wrapping the nucleotide with the given offset
     */
    public Optional<Nucleotide> getNucleotideWithOffset(int offset) {
        try {
            Chain chain = getParentChain();
            int index = chain.getGroups().indexOf(this);
            return Optional.of((Nucleotide) chain.getGroups().get(index + offset));
        } catch (Exception e) {
            return Optional.empty();
        }
    }

    /**
     * Access to the next nucleotide in the chain.
     * <code>getPreviousNucleotide() is equivalent to
     * </code><blockquote><pre>
     * getNucleotideWithOffset(1)
     * </pre></blockquote>
     * @see #getNucleotideWithOffset(int)
     * @return an optional wrapping the next nucleotide
     */
    public Optional<Nucleotide> getNextNucleotide() {
        return getNucleotideWithOffset(1);
    }
}
