package de.bioforscher.jstructure.graph;

import de.bioforscher.jstructure.feature.interaction.HydrogenBond;
import de.bioforscher.jstructure.feature.interaction.HydrophobicInteraction;
import de.bioforscher.jstructure.feature.interaction.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.feature.sse.SecondaryStructureType;
import de.bioforscher.jstructure.mathematics.Pair;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;

import java.util.ArrayList;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.function.BiPredicate;
import java.util.stream.Collectors;

public class ResidueGraph extends SimpleGraph<AminoAcid, DefaultEdge> {
    private final List<AminoAcid> aminoAcids;
    private final List<Pair<AminoAcid, AminoAcid>> contacts;
    private final List<Pair<AminoAcid, AminoAcid>> localContacts;
    private final List<Pair<AminoAcid, AminoAcid>> nonLocalContacts;
    private final String sequence;
    private final String secondaryStructureElements;

    public enum InteractionScheme {
        CALPHA8((aminoAcid1, aminoAcid2) -> aminoAcid1.getCa().calculate()
                .distanceFast(aminoAcid2.getCa()) < 8 * 8),
        SALENTIN2015((aminoAcid1, aminoAcid2) -> aminoAcid1.getParentChain().getFeature(PLIPInteractionContainer.class)
                .areInContact(aminoAcid1, aminoAcid2)),
        SALENTIN2015_HYDROGEN_BONDS(((aminoAcid1, aminoAcid2) -> aminoAcid1.getParentChain().getFeature(PLIPInteractionContainer.class)
                .areInContactByInteractionType(aminoAcid1, aminoAcid2, HydrogenBond.class))),
        SALENTIN2015_HYDROPHOBIC_INTERACTION(((aminoAcid1, aminoAcid2) -> aminoAcid1.getParentChain().getFeature(PLIPInteractionContainer.class)
                .areInContactByInteractionType(aminoAcid1, aminoAcid2, HydrophobicInteraction.class)));

        private final BiPredicate<AminoAcid, AminoAcid> criterion;

        InteractionScheme(BiPredicate<AminoAcid, AminoAcid> criterion) {
            this.criterion = criterion;
        }

        public boolean areInContact(AminoAcid aminoAcid1, AminoAcid aminoAcid2) {
            return criterion.test(aminoAcid1, aminoAcid2);
        }
    }

    public static ResidueGraph createResidueGraph(Chain chain, InteractionScheme interactionScheme) {
        List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());
        List<Pair<AminoAcid, AminoAcid>> contacts = new ArrayList<>();

        for(int i = 0; i < aminoAcids.size() - 1; i++) {
            AminoAcid aminoAcid1 = aminoAcids.get(i);
            for(int j = i + 1; j < aminoAcids.size(); j++) {
                AminoAcid aminoAcid2 = aminoAcids.get(j);
                if(j == i + 1) {
                    contacts.add(new Pair<>(aminoAcid1, aminoAcid2));
                    continue;
                }

                if(interactionScheme.areInContact(aminoAcid1, aminoAcid2)) {
                    contacts.add(new Pair<>(aminoAcid1, aminoAcid2));
                }
            }
        }

        return new ResidueGraph(aminoAcids, contacts);
    }

    public ResidueGraph(List<AminoAcid> aminoAcids, List<Pair<AminoAcid, AminoAcid>> contacts) {
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
        this.sequence = aminoAcids.stream()
                .map(AminoAcid::getOneLetterCode)
                .collect(Collectors.joining());
        this.secondaryStructureElements = aminoAcids.stream()
                .map(aminoAcid -> aminoAcid.getFeature(GenericSecondaryStructure.class))
                .map(GenericSecondaryStructure::getSecondaryStructure)
                .map(SecondaryStructureType::getReducedRepresentation)
                .collect(Collectors.joining());
    }

    public String getSequence() {
        return sequence;
    }

    public String getSecondaryStructureElements() {
        return secondaryStructureElements;
    }

    public String getCaspRRRepresentation() {
        // i  j  d1  d2  p
        // 24 33 0 8 12.279515
        // i and j: residue numbers
        return contacts.stream()
                .map(edge -> (edge.getLeft().getResidueIndex() + 1) + " " +
                        (edge.getRight().getResidueIndex() + 1) + " 0 8.0 1")
                .collect(Collectors.joining(System.lineSeparator(),
                        sequence + System.lineSeparator(),
                        ""));
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
