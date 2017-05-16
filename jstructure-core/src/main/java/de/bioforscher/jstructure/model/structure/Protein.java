package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.model.feature.AbstractFeatureContainer;
import de.bioforscher.jstructure.model.structure.container.ChainContainer;
import de.bioforscher.jstructure.model.structure.selection.Selection;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * The root element of protein structures. A really simplified de.bioforscher.explorer.helices.model, reduced to our particularized needs. Consists of
 * {@link Chain} objects.<br />
 * Created by S on 27.09.2016.
 */
public class Protein extends AbstractFeatureContainer implements ChainContainer {
    /**
     * reference to an undefined protein - this is used by chains without explicit parent reference
     */
    static final Protein UNKNOWN_PROTEIN = new Protein();
    static {
        UNKNOWN_PROTEIN.setName("UNKNOWN PROTEIN");
    }
    /**
     * The <tt>PDB</tt> id of this protein. If none is known (e.g. because this is a modeled structure, the filename is returned)
     */
    private String name;
    /**
     * The <tt>PDB</tt> description of this protein.
     */
    private String title;
    private List<Chain> chains;
    private String identifier;

    public Protein() {
        this.chains = new ArrayList<>();
    }

    public Protein(Protein protein) {
        this.name = protein.name;
        this.title = protein.title;
        this.chains = protein.chains()
                .map(Chain::new)
                .collect(Collectors.toList());
        this.chains.forEach(chain -> chain.setParentProtein(this));
        this.identifier = protein.identifier;
    }

    public Protein(List<Chain> chains) {
        this.chains = chains;
    }

    public Selection.ChainSelection select() {
        return Selection.on(this);
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
     * The number of groups associated to this protein.
     * @return the number of amino acids, hetatms and nucleotids
     */
    public int getSize() {
        return (int) groups().count();
    }

    /**
     * The name of this structure as <tt>PDB</tt> id (or the parsed file's name as fallback).
     * @return the name of this protein
     */
    public String getName() {
        return name;
    }

    /**
     * Assign a name (i.e. most of the time a <tt>PDB</tt> id) to this protein.
     * @param name the <tt>PDB</tt> id found in the file which was parsed to create this protein, otherwise the filename,
     *             otherwise probably <code>null</code>
     */
    public void setName(String name) {
        this.name = name;
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

    /**
     * Registers a child. This object will assign a reference to itself to the getChain.
     * @param chain the getChain to process
     */
    public void addChain(Chain chain) {
        getChains().add(chain);
        chain.setParentProtein(this);
    }

    @Override
    public String getIdentifier() {
        return identifier == null ? name : identifier;
    }

    @Override
    public String toString() {
        return getClass().getSimpleName() + " identifier='" + getIdentifier() + "' chains='"  + getChains().size() + "'";
    }

    public void setIdentifier(String identifier) {
        this.identifier = identifier;
    }
}