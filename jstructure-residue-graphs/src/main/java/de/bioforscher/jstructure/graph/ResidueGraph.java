package de.bioforscher.jstructure.graph;

import de.bioforscher.jstructure.graph.contact.definition.ContactDefinition;
import de.bioforscher.jstructure.graph.contact.definition.ContactDefinitionFactory;
import de.bioforscher.jstructure.mathematics.Pair;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;

import java.util.ArrayList;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.stream.Collectors;

/**
 * Represents a protein as residue graph. Covalent bonds and local contacts are considered too.
 */
public class ResidueGraph extends SimpleGraph<AminoAcid, DefaultEdge> {
    protected final List<AminoAcid> aminoAcids;
    protected final List<Pair<AminoAcid, AminoAcid>> contacts;
    protected final List<Pair<AminoAcid, AminoAcid>> localContacts;
    protected final List<Pair<AminoAcid, AminoAcid>> longRangeContacts;
    protected final ContactDefinition contactDefinition;

    public static ResidueGraph createResidueGraph(Chain chain, ContactDefinition contactDefinition) {
        List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());
        List<Pair<AminoAcid, AminoAcid>> contacts = createContactList(aminoAcids, contactDefinition);
        return new ResidueGraph(aminoAcids, contacts, contactDefinition);
    }

    protected static List<Pair<AminoAcid, AminoAcid>> createContactList(List<AminoAcid> aminoAcids, ContactDefinition contactDefinition) {
        List<Pair<AminoAcid, AminoAcid>> contacts = new ArrayList<>();

        for(int i = 0; i < aminoAcids.size() - 1; i++) {
            AminoAcid aminoAcid1 = aminoAcids.get(i);
            for(int j = i + 1; j < aminoAcids.size(); j++) {
                AminoAcid aminoAcid2 = aminoAcids.get(j);
                if(j == i + 1) {
                    contacts.add(new Pair<>(aminoAcid1, aminoAcid2));
                    continue;
                }

                if(contactDefinition.areInContact(aminoAcid1, aminoAcid2)) {
                    contacts.add(new Pair<>(aminoAcid1, aminoAcid2));
                }
            }
        }

        return contacts;
    }

    public static ResidueGraph createPlipResidueGraph(Chain chain) {
        return createResidueGraph(chain, ContactDefinitionFactory.createPlipContactDefinition());
    }

    public static ResidueGraph createHydrogenBondPlipResidueGraph(Chain chain) {
        return createResidueGraph(chain, ContactDefinitionFactory.createHydrogenBondContactDefinition());
    }

    public static ResidueGraph createHydrophobicInteractionResidueGraph(Chain chain) {
        return createResidueGraph(chain, ContactDefinitionFactory.createHydrophobicInteractionContactDefinition());
    }

    public static ResidueGraph createDistanceResidueGraph(Chain chain) {
        return createResidueGraph(chain, ContactDefinitionFactory.createAlphaCarbonContactDefinition(8.0));
    }

    ResidueGraph(List<AminoAcid> aminoAcids, List<Pair<AminoAcid, AminoAcid>> contacts, ContactDefinition contactDefinition) {
        super(DefaultEdge.class);
        this.aminoAcids = aminoAcids;
        aminoAcids.forEach(this::addVertex);
        this.contacts = contacts;
        contacts.forEach(pair -> addEdge(pair.getLeft(), pair.getRight()));
        this.localContacts = contacts.stream()
                .filter(ResidueGraph::isLocalContact)
                .collect(Collectors.toList());
        this.longRangeContacts = contacts.stream()
                .filter(ResidueGraph::isLongRangeContact)
                .collect(Collectors.toList());
        this.contactDefinition = contactDefinition;
    }

    public ContactDefinition getContactDefinition() {
        return contactDefinition;
    }

    public String getConfoldRRType() {
        return contactDefinition.getConfoldRRType();
    }

    protected static boolean isLongRangeContact(Pair<AminoAcid, AminoAcid> pair) {
        return Math.abs(pair.getLeft().getResidueIndex() - pair.getRight().getResidueIndex()) > 5;
    }

    private static boolean isLocalContact(Pair<AminoAcid, AminoAcid> pair) {
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

    public List<Pair<AminoAcid, AminoAcid>> getLongRangeContacts() {
        return longRangeContacts;
    }

    public List<AminoAcid> getNonLocalContactsOf(AminoAcid aminoAcid) {
        return filterContactList(longRangeContacts, aminoAcid);
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
