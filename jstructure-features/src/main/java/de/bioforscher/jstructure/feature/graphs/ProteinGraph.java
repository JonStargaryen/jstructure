package de.bioforscher.jstructure.feature.graphs;

import de.bioforscher.jstructure.mathematics.Pair;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;

import java.util.List;
import java.util.NoSuchElementException;
import java.util.stream.Collectors;

public class ProteinGraph extends SimpleGraph<AminoAcid, DefaultEdge> {
    private final List<AminoAcid> aminoAcids;
    private final List<Pair<AminoAcid, AminoAcid>> contacts;
    private final List<Pair<AminoAcid, AminoAcid>> localContacts;
    private final List<Pair<AminoAcid, AminoAcid>> nonLocalContacts;

    public ProteinGraph(List<AminoAcid> aminoAcids, List<Pair<AminoAcid, AminoAcid>> contacts) {
        super(DefaultEdge.class);
        this.aminoAcids = aminoAcids;
        aminoAcids.forEach(this::addVertex);
        this.contacts = contacts;
        contacts.forEach(pair -> addEdge(pair.getLeft(), pair.getRight()));
        this.localContacts = contacts.stream()
                .filter(this::isLocalContact)
                .collect(Collectors.toList());
        this.nonLocalContacts = contacts.stream()
                .filter(this::isNonLocalContact)
                .collect(Collectors.toList());
    }

    private boolean isNonLocalContact(Pair<AminoAcid, AminoAcid> pair) {
        return Math.abs(pair.getLeft().getResidueIndex() - pair.getRight().getResidueIndex()) > 5;
    }

    private boolean isLocalContact(Pair<AminoAcid, AminoAcid> pair) {
        return Math.abs(pair.getLeft().getResidueIndex() - pair.getRight().getResidueIndex()) < 6;
    }

    public List<AminoAcid> getAminoAcids() {
        return aminoAcids;
    }

    public List<Pair<AminoAcid, AminoAcid>> getContacts() {
        return contacts;
    }

    public List<AminoAcid> getContactsOf(AminoAcid aminoAcid) {
        return filterContactList(contacts, aminoAcid);
    }

    public List<AminoAcid> getContactsOf(ResidueIdentifier residueIdentifier) {
        return getContactsOf(resolve(residueIdentifier));
    }

    private List<AminoAcid> filterContactList(List<Pair<AminoAcid, AminoAcid>> list, AminoAcid aminoAcid) {
        return list.stream()
                .filter(pair -> filterPair(pair, aminoAcid))
                .map(pair -> resolvePair(pair, aminoAcid))
                .collect(Collectors.toList());
    }

    private boolean filterPair(Pair<AminoAcid, AminoAcid> pair, AminoAcid aminoAcid) {
        return pair.contains(aminoAcid);
    }

    private AminoAcid resolvePair(Pair<AminoAcid, AminoAcid> pair, AminoAcid aminoAcid) {
        return pair.getLeft().equals(aminoAcid) ? pair.getRight() : pair.getLeft();
    }

    public List<Pair<AminoAcid, AminoAcid>> getLocalContacts() {
        return localContacts;
    }

    public List<AminoAcid> getLocalContactsOf(AminoAcid aminoAcid) {
        return filterContactList(localContacts, aminoAcid);
    }

    public List<AminoAcid> getLocalContactsOf(ResidueIdentifier residueIdentifier) {
        return getLocalContactsOf(resolve(residueIdentifier));
    }

    public List<Pair<AminoAcid, AminoAcid>> getNonLocalContacts() {
        return nonLocalContacts;
    }

    public List<AminoAcid> getNonLocalContactsOf(AminoAcid aminoAcid) {
        return filterContactList(nonLocalContacts, aminoAcid);
    }

    public List<AminoAcid> getNonLocalContactsOf(ResidueIdentifier residueIdentifier) {
        return getNonLocalContactsOf(resolve(residueIdentifier));
    }

    private AminoAcid resolve(ResidueIdentifier residueIdentifier) {
        return aminoAcids.stream()
                .filter(aminoAcid -> aminoAcid.getResidueIdentifier().equals(residueIdentifier))
                .findFirst()
                .orElseThrow(() ->  new NoSuchElementException("did not find amino acid with id " + residueIdentifier));
    }
}
