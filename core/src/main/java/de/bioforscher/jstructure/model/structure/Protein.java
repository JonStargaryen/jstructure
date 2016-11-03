package de.bioforscher.jstructure.model.structure;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * The root element of protein structures. A really simplified model, reduced to our particularized needs. Consists of
 * {@link Chain} objects.<br />
 * Created by S on 27.09.2016.
 * TODO: at some point a wide range of fields should probably become final
 */
public class Protein implements ChainContainer, AtomRecordWriter {
    /**
     * The chain container.
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

    @Override
    public void clear() {
        chains().forEach(Chain::clear);
    }

    @Override
    public void clearPseudoAtoms() {
        chains.forEach(Chain::clearPseudoAtoms);
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
     * Registers a child. This object will assign a reference to itself to the chain.
     * @param chain the chain to process
     */
    public void addChain(Chain chain) {
        chains.add(chain);
        chain.setParentProtein(this);
    }

    @Override
    public Stream<Chain> chains() {
        return chains.stream();
    }

    @Override
    public Stream<Group> groups() {
        return chains().flatMap(Chain::groups);
    }

    @Override
    public Stream<Residue> residues() {
        return chains().flatMap(Chain::residues);
    }

    @Override
    public Stream<Atom> atoms() {
        return groups().flatMap(Group::atoms);
    }

    @Override
    public String toString() {
        return getClass().getSimpleName() + " name='" + name + "'";
    }

    @Override
    public String composePDBRecord() {
        return chains.stream().map(Chain::composePDBRecord).collect(Collectors.joining(System.lineSeparator()));
    }

    @Override
    public Map<Enum, Object> getFeatureMap() {
        return featureMap;
    }
}