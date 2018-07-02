package de.bioforscher.jstructure.model.structure.nucleotide;

import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.*;

import java.lang.reflect.Field;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Representation of a nucleotide.
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
    private static final Set<String> ASSIGNABLE_ATOM_NAMES = Stream.of("op3",
            "p",
            "op1",
            "op2",
            "o5prime",
            "c5prime",
            "c4prime",
            "o4prime",
            "c3prime",
            "o3prime",
            "c2prime",
            "c1prime",
            "n1",
            "c2",
            "n3",
            "c4",
            "c5",
            "c6")
            .collect(Collectors.toSet());

    //TODO grouping of atom names & family enum

    Nucleotide(Nucleotide nucleotide, boolean deep) {
        super(nucleotide, deep);
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
                return new UnknownNucleotide(pdbName,
                        residueIdentifier,
                        ligand);
            } else {
                try {
                    return representingClass.getConstructor(ResidueIdentifier.class, boolean.class)
                            .newInstance(residueIdentifier, ligand);
                } catch (Exception e) {
                    throw new RuntimeException("creation of Nucleotide instance failed", e);
                }
            }
        }

        public static Nucleotide.Family resolveOneLetterCode(char oneLetterCode) {
            return resolveOneLetterCode(String.valueOf(oneLetterCode));
        }

        public static Nucleotide.Family resolveOneLetterCode(String oneLetterCode) {
            return Stream.of(Nucleotide.Family.values())
                    .filter(nucleotide -> oneLetterCode.equalsIgnoreCase(nucleotide.getOneLetterCode()))
                    .findFirst()
                    .orElse(Family.UNKNOWN_NUCLEOTIDE);
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

        public String getOneLetterCode() {
            return groupPrototype.getOneLetterCode().get();
        }

        @Override
        public String getThreeLetterCode() {
            return groupPrototype.getThreeLetterCode();
        }
    }

    public Optional<Atom> getOp3Optional() {
        return Optional.ofNullable(op3);
    }

    public Atom getOp3() {
        return getOp3Optional().orElseThrow(() -> createNoAtomException("OP3"));
    }

    public Optional<Atom> getPOptional() {
        return Optional.ofNullable(p);
    }

    public Atom getP() {
        return getPOptional().orElseThrow(() -> createNoAtomException("P"));
    }

    public Optional<Atom> getOp1Optional() {
        return Optional.ofNullable(op1);
    }

    public Atom getOp1() {
        return getOp1Optional().orElseThrow(() -> createNoAtomException("OP1"));
    }

    public Optional<Atom> getOp2Optional() {
        return Optional.ofNullable(op2);
    }

    public Atom getOp2() {
        return getOp2Optional().orElseThrow(() -> createNoAtomException("OP2"));
    }

    public Optional<Atom> getO5primeOptional() {
        return Optional.ofNullable(o5prime);
    }

    public Atom getO5prime() {
        return getO5primeOptional().orElseThrow(() -> createNoAtomException("O5'"));
    }

    public Optional<Atom> getC5primeOptional() {
        return Optional.ofNullable(c5prime);
    }

    public Atom getC5prime() {
        return getC5primeOptional().orElseThrow(() -> createNoAtomException("C5'"));
    }

    public Optional<Atom> getC4primeOptional() {
        return Optional.ofNullable(c4prime);
    }

    public Atom getC4prime() {
        return getC4primeOptional().orElseThrow(() -> createNoAtomException("C4'"));
    }

    public Optional<Atom> getO4primeOptional() {
        return Optional.ofNullable(o4prime);
    }

    public Atom getO4prime() {
        return getO4primeOptional().orElseThrow(() -> createNoAtomException("O4'"));
    }

    public Optional<Atom> getC3primeOptional() {
        return Optional.ofNullable(c3prime);
    }

    public Atom getC3prime() {
        return getC3primeOptional().orElseThrow(() -> createNoAtomException("C3'"));
    }

    public Optional<Atom> getO3primeOptional() {
        return Optional.ofNullable(o3prime);
    }

    public Atom getO3prime() {
        return getO3primeOptional().orElseThrow(() -> createNoAtomException("O3'"));
    }

    public Optional<Atom> getC2primeOptional() {
        return Optional.ofNullable(c2prime);
    }

    public Atom getC2prime() {
        return getC2primeOptional().orElseThrow(() -> createNoAtomException("C2'"));
    }

    public Optional<Atom> getC1primeOptional() {
        return Optional.ofNullable(c1prime);
    }

    public Atom getC1prime() {
        return getC1primeOptional().orElseThrow(() -> createNoAtomException("C1'"));
    }

    public Optional<Atom> getN1Optional() {
        return Optional.ofNullable(n1);
    }

    public Atom getN1() {
        return getN1Optional().orElseThrow(() -> createNoAtomException("N1"));
    }

    public Optional<Atom> getC2Optional() {
        return Optional.ofNullable(c2);
    }

    public Atom getC2() {
        return getC2Optional().orElseThrow(() -> createNoAtomException("C2"));
    }

    public Optional<Atom> getN3Optional() {
        return Optional.ofNullable(n3);
    }

    public Atom getN3() {
        return getN3Optional().orElseThrow(() -> createNoAtomException("N3"));
    }

    public Optional<Atom> getC4Optional() {
        return Optional.ofNullable(c4);
    }

    public Atom getC4() {
        return getC4Optional().orElseThrow(() -> createNoAtomException("C4"));
    }

    public Optional<Atom> getC5Optional() {
        return Optional.ofNullable(c5);
    }

    public Atom getC5() {
        return getC5Optional().orElseThrow(() -> createNoAtomException("C5"));
    }

    public Optional<Atom> getC6Optional() {
        return Optional.ofNullable(c6);
    }

    public Atom getC6() {
        return getC6Optional().orElseThrow(() -> createNoAtomException("C6"));
    }

    public String getOneLetterCode() {
        return getGroupPrototype().getOneLetterCode().orElse("?");
    }

    @Override
    protected void addAtomInternal(Atom atom) {
        // assign atoms by reflection
        String atomName = atom.getName();
        String fieldName = atomName.toLowerCase()
                // remove additional characters
                .replace("\"", "")
                // replace ' by prime
                .replace("'", "prime");

        // some field names will never be presented: e.g. hydrogen atom
        if(fieldName.contains("h") && !fieldName.equals("h")) {
            return;
        }

        try {
            // find field in Nucleotide or child classes
            Field field = ASSIGNABLE_ATOM_NAMES.contains(fieldName) ?
                    this.getClass().getSuperclass().getDeclaredField(fieldName) :
                    this.getClass().getDeclaredField(fieldName);
            field.setAccessible(true);
            Object fieldContent = field.get(this);

            // protect already set fields
            if (fieldContent != null) {
                return;
            }

            field.set(this, atom);
        } catch (NoSuchFieldException | IllegalAccessException e) {
            if(this.isStandardNucleotide() && !this.getClass().equals(UnknownNucleotide.class)) {
                logger.warn("missing field for atom {} in class {}",
                        atom.getName(),
                        this.getClass().getSimpleName(),
                        e);
            }
        }
    }

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
