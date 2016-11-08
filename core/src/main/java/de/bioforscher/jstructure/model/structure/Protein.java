package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.model.structure.container.ChainContainer;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * The root element of protein structures. A really simplified model, reduced to our particularized needs. Consists of
 * {@link Chain} objects.<br />
 * Created by S on 27.09.2016.
 */
public class Protein implements ChainContainer, AtomRecordWriter {
    /**
     * The getChain container.
     */
    private List<Chain> chains;
    /**
     * The <tt>PDB</tt> id of this protein. If none is known (e.g. because this is a modeled structure, the filename is returned)
     */
    private String name;
    /**
     * The <tt>PDB</tt> description of this protein.
     */
    private String title;
    private Map<Enum, Object> featureMap;

    /**
     * The default constructor.
     */
    public Protein() {
        this.chains = new ArrayList<>();
        this.featureMap = new HashMap<>();
    }

    public List<Chain> getChains() {
        return chains;
    }

    public List<Group> getGroups() {
        return chains()
                .flatMap(Chain::groups)
                .collect(Collectors.toList());
    }

    public List<Atom> getAtoms() {
        return getGroups().stream()
                .flatMap(Group::atoms)
                .collect(Collectors.toList());
    }

    /**
     * The number of groups associated to this protein.
     * @return the number of amino acids, hetatms and nucleotids
     */
    public int getSize() {
        return (int) residues().count();
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
        chains.add(chain);
        chain.setParentProtein(this);
    }

    @Override
    public String toString() {
        return getClass().getSimpleName() + " name='" + name + "'";
    }

    @Override
    public Map<Enum, Object> getFeatureMap() {
        return featureMap;
    }
}