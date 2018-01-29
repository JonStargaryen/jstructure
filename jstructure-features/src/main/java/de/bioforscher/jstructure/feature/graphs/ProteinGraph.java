package de.bioforscher.jstructure.feature.graphs;

import de.bioforscher.jstructure.mathematics.Pair;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;

import java.util.List;
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
        return Math.abs(pair.getLeft().getResidueIndex() - pair.getRight().getResidueIndex()) > 6;
    }

    private boolean isLocalContact(Pair<AminoAcid, AminoAcid> pair) {
        return Math.abs(pair.getLeft().getResidueIndex() - pair.getRight().getResidueIndex()) <= 6;
    }

    public List<AminoAcid> getAminoAcids() {
        return aminoAcids;
    }

    public List<Pair<AminoAcid, AminoAcid>> getContacts() {
        return contacts;
    }

    public List<Pair<AminoAcid, AminoAcid>> getLocalContacts() {
        return localContacts;
    }

    public List<Pair<AminoAcid, AminoAcid>> getNonLocalContacts() {
        return nonLocalContacts;
    }
}
