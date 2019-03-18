package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.mathematics.LinearAlgebra;
import de.bioforscher.jstructure.model.feature.AbstractFeatureable;
import de.bioforscher.jstructure.model.identifier.ProteinIdentifier;
import de.bioforscher.jstructure.model.structure.container.ChainContainer;
import de.bioforscher.jstructure.model.structure.selection.Selection;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * The root element of protein structures. A really simplified de.bioforscher.explorer.helices.model, reduced to our particularized needs. Consists of
 * {@link Chain} objects.<br />
 * Created by S on 27.09.2016.
 */
public class Structure extends AbstractFeatureable implements ChainContainer/*, FeatureContainerRoot*/ {
    /**
     * reference to an undefined protein - this is used by chains without explicit parent reference
     */
    public static final Structure UNKNOWN_STRUCTURE = new Structure(ProteinIdentifier.UNKNOWN_PROTEIN_IDENTIFIER);
    private ProteinIdentifier proteinIdentifier;
    /**
     * The <tt>PDB</tt> description of this protein.
     */
    private String title;
    /**
     * The shorter String in the <tt>PDB</tt> header.
     * HEADER    COMPLEX (OXIDOREDUCTASE/ANTIBODY)       08-AUG-97   1AR1
     */
    private String classification;
    private LocalDate depositionDate;
    private List<Chain> chains;
    private String identifier;

    public Structure(ProteinIdentifier proteinIdentifier) {
        this.proteinIdentifier = proteinIdentifier;
        this.classification = proteinIdentifier.getAdditionalName().isEmpty() ?
                ProteinIdentifier.UNKNOWN_PROTEIN_IDENTIFIER.getAdditionalName() : proteinIdentifier.getAdditionalName();
        this.depositionDate = LocalDate.now();
        this.title = classification;
        this.chains = new ArrayList<>();
    }

    Structure(Structure structure, boolean deepCopy) {
        this.proteinIdentifier = structure.proteinIdentifier;
        this.classification = structure.classification;
        this.depositionDate = structure.depositionDate;
        this.title = structure.title;
        this.identifier = structure.identifier;
        if(deepCopy) {
            this.chains = structure.chains()
                    .map(Chain::createDeepCopy)
                    .collect(Collectors.toList());
            this.chains.forEach(chain -> chain.setParentStructure(this));
        } else {
            this.chains = new ArrayList<>();
        }
    }

    public Selection.ChainSelection select() {
        return Selection.on(this);
    }

    public LinearAlgebra.AtomContainerLinearAlgebra calculate() {
        return LinearAlgebra.on(this);
    }

    public Chain getFirstChain() {
        return chains.get(0);
    }

    @Override
    public List<Chain> getChains() {
        return chains;
    }

    @Override
    public List<Group> getGroups() {
        return chains().flatMap(Chain::groups).collect(Collectors.toList());
    }

    @Override
    public List<Atom> getAtoms() {
        return chains().flatMap(Chain::atoms).collect(Collectors.toList());
    }

    /**
     * The name of this structure as <tt>PDB</tt> id (or the parsed file's name as fallback).
     * @return the id of this protein
     */
    public ProteinIdentifier getProteinIdentifier() {
        return proteinIdentifier;
    }

    /**
     * Assign a name (i.e. most of the time a <tt>PDB</tt> id) to this protein.
     * @param pdbId the <tt>PDB</tt> id found in the file which was parsed to create this protein, otherwise the filename,
     *             otherwise probably <code>null</code>
     */
    void setProteinIdentifier(ProteinIdentifier pdbId) {
        this.proteinIdentifier = pdbId;
    }

    @Override
    public Structure getParentStructure() {
        return this;
    }

    /**
     * The description of this protein, e.g. found in <tt>PDB</tt> headers.
     * @return the protein title as String
     */
    public String getTitle() {
        return title;
    }

    /**
     * Assign a description to this protein.
     * @param title a more or less detailed description of this protein
     */
    public void setTitle(String title) {
        this.title = title;
    }

    public String getClassification() {
        return classification;
    }

    public void setClassification(String classification) {
        this.classification = classification;
    }

    public LocalDate getDepositionDate() {
        return depositionDate;
    }

    public void setDepositionDate(LocalDate depositionDate) {
        this.depositionDate = depositionDate;
    }

    public String getHeader() {
        // formatting by BioJava code as usual
        StringBuilder stringBuilder = new StringBuilder();

        // too long classification strings could get introduced
        String classificationString = classification.length() < 40 ? classification : classification.substring(0, 40);
        stringBuilder.append("HEADER    ");
        stringBuilder.append(classificationString);
//        stringBuilder.append(" ");

        // fill up the white space to the right column
        int l =  classificationString.length() + 10;
        while (l < 50) {
            stringBuilder.append(" ");
            l++;
        }

        stringBuilder.append(StandardFormat.format(depositionDate));
        stringBuilder.append("   ");

        String pdbId;
        if(getProteinIdentifier() != null && getProteinIdentifier().getPdbId() != null) {
            pdbId = getProteinIdentifier().getPdbId().toUpperCase();
        } else {
            pdbId = "1XXX";
        }
        stringBuilder.append(pdbId);

        // final padding
        while (stringBuilder.length() < 80) {
            stringBuilder.append(" ");
        }

        stringBuilder.append(System.lineSeparator());
        return stringBuilder.toString();
    }

    /**
     * Registers a child. This object will assign a reference to itself to the getChain.
     * @param chain the getChain to processUniProtId
     */
    public void addChain(Chain chain) {
        getChains().add(chain);
        chain.setParentStructure(this);
    }

    @Override
    public String toString() {
        return getClass().getSimpleName() + " " + getIdentifier() + " chains="  + getChains().size();
    }

    @Override
    public String getIdentifier() {
        return identifier == null ? proteinIdentifier.getFullName() : identifier;
    }

    @Override
    public void setIdentifier(String identifier) {
        this.identifier = identifier;
    }

    @Override
    public Structure createDeepCopy() {
        return new Structure(this, true);
    }

    @Override
    public Structure createShallowCopy() {
        return new Structure(this, false);
    }
}