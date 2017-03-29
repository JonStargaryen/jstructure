package de.bioforscher.jstructure.model.structure.family;

import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * The container of parsed group information.
 * Created by S on 06.01.2017.
 */
public class GroupInformation {
    public static final GroupInformation UNKNOWN_GROUP = builder()
            .name("UNKNOWN GROUP")
            .oneLetterCode("X")
            .threeLetterCode("UNK")
            .type("NON-POLYMER")
            .build();

    public static final GroupInformation UNKNOWN_AMINO_ACID = builder()
            .name("UNKNOWN AMINO ACID")
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

    GroupInformation() {

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
            // fix 17/03/29: names on multiple lines are not parsed - they contain lots of white-spaces and colons
            String document = lines.collect(Collectors.joining(System.lineSeparator()));
            if(document.contains("_chem_comp.name")) {
                name(splitCifDocument(document, "_chem_comp.name"));
            }
            if(document.contains("_chem_comp.type")) {
                type(splitCifDocument(document, "_chem_comp.type"));
            }
            if(document.contains("_chem_comp.one_letter_code")) {
                oneLetterCode(splitCifDocument(document, "_chem_comp.one_letter_code"));
            }
            if(document.contains("_chem_comp.three_letter_code")) {
                threeLetterCode(splitCifDocument(document, "_chem_comp.three_letter_code"));
            }
            if(document.contains("_chem_comp.mon_nstd_parent_comp_id")) {
                parentCompound(splitCifDocument(document, "_chem_comp.mon_nstd_parent_comp_id"));
            }
            return this.build();
        }

        // fix 17/03/29: adds long-entry support
        private String splitCifDocument(String document, String tag) {
            return document.split(tag)[1].split("_")[0].replaceAll("[\\r\\n]", "").replace(";", "").trim();
        }

         @Deprecated
        private String splitCifLine(String line) {
            return line.substring(49).trim();
        }
    }
}
