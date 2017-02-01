package de.bioforscher.jstructure.model.structure.family;

import java.util.stream.Stream;

/**
 * The container of parsed group information.
 * Created by S on 06.01.2017.
 */
public class GroupInformation {
    public static final GroupInformation UNKNOWN_AMINO_ACID = builder()
            .name("UNKNOWN GROUP")
            .oneLetterCode("X")
            .threeLetterCode("UNK")
            .type("PEPTIDE LINKING")
            .build();

    /**
     * The full name of this group.
     */
    private String name;
    /**
     * The one-letter-code of this group.
     */
    private String oneLetterCode;
    /**
     * The three-letter-code (or PDB name) of this group.
     */
    private String threeLetterCode;
    /**
     * The raw type of this group (i.e. amino acid : PEPTIDE LINKING, nucleotide : RNA LINKING or ligand : NON-POLYMER).
     */
    private String type;
    /**
     * The parent of e.g. modified amino acids.
     */
    private String parentCompound;

    private GroupInformation(Builder builder) {
        this.name = builder.name;
        this.oneLetterCode = builder.oneLetterCode;
        this.threeLetterCode = builder.threeLetterCode;
        this.type = builder.type;
        this.parentCompound = builder.parentCompound;
    }

    public String getName() {
        return name;
    }

    public String getOneLetterCode() {
        return oneLetterCode;
    }

    public String getThreeLetterCode() {
        return threeLetterCode;
    }

    public String getType() {
        return type;
    }

    public String getParentCompound() {
        return parentCompound;
    }

    public AminoAcidFamily getAminoAcidFamily() {
        return AminoAcidFamily.valueOfIgnoreCase(threeLetterCode).orElse(AminoAcidFamily.UNKNOWN);
    }

    public NucleotideFamily getNucleotideFamily() {
        return NucleotideFamily.valueOfIgnoreCase(threeLetterCode).orElse(NucleotideFamily.UNKNOWN);
    }

    public boolean isAminoAcid() {
        //TODO this could fail for single amino acids in a structure which are actually ligands
        return type.contains("PEPTIDE LINKING") && (!getAminoAcidFamily().equals(AminoAcidFamily.UNKNOWN) || isModified());
    }

    public boolean isNucleotide() {
        //TODO this could fail for single nucleotides in a structure which are actually ligands
        return type.contains("NA LINKING");
    }

    public boolean isLigand() {
        return !isAminoAcid() && !isNucleotide();
    }

    public boolean isModified() {
        return !parentCompound.equals(Builder.UNSET_PARENT_COMPOUND);
    }

    @Override
    public String toString() {
        return getClass().getSimpleName() + " name='" + name + "' one-letter-code='" + oneLetterCode + "' three-letter-code='" +
                threeLetterCode + "' type='" + type + "'";
    }

    public static Builder builder() {
        return new Builder();
    }

    /**
     * The builder based on CIF files.
     */
    public static class Builder {
        String name;
        String oneLetterCode;
        String threeLetterCode;
        String type;
        String parentCompound = UNSET_PARENT_COMPOUND;
        static final String UNSET_PARENT_COMPOUND = "?";

        public GroupInformation build() {
            return new GroupInformation(this);
        }

        public Builder name(String name) {
            this.name = name;
            return this;
        }

        public Builder oneLetterCode(String oneLetterCode) {
            this.oneLetterCode = oneLetterCode;
            return this;
        }

        public Builder threeLetterCode(String threeLetterCode) {
            this.threeLetterCode = threeLetterCode;
            return this;
        }

        public Builder type(String type) {
            this.type = type;
            return this;
        }

        public Builder parentCompound(String parentCompound) {
            this.parentCompound = parentCompound;
            return this;
        }

        public GroupInformation parseLines(Stream<String> lines) {
            lines.forEach(line -> {
                if(line.startsWith("_chem_comp.name")) {
                    name(splitCifLine(line));
                }
                if(line.startsWith("_chem_comp.type")) {
                    type(splitCifLine(line));
                }
                if(line.startsWith("_chem_comp.one_letter_code")) {
                    oneLetterCode(splitCifLine(line));
                }
                if(line.startsWith("_chem_comp.three_letter_code")) {
                    threeLetterCode(splitCifLine(line));
                }
                if(line.startsWith("_chem_comp.mon_nstd_parent_comp_id")) {
                    parentCompound(splitCifLine(line));
                }
            });
            return this.build();
        }

        private String splitCifLine(String line) {
            return line.substring(49).trim();
        }
    }
}
