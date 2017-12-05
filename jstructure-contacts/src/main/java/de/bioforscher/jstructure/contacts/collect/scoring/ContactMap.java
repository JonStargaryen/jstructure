package de.bioforscher.jstructure.contacts.collect.scoring;

import de.bioforscher.jstructure.contacts.ContactsConstants;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.mathematics.Pair;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

public class ContactMap {
    private final List<AminoAcid> aminoAcids;
    private final List<Contact> contacts;

    private ContactMap(List<AminoAcid> aminoAcids, List<Contact> contacts) {
        this.aminoAcids = aminoAcids;
        this.contacts = contacts;
    }

    public static ContactMap createPlipContactMap(Structure structure) {
        List<AminoAcid> aminoAcids = structure.aminoAcids().collect(Collectors.toList());
        List<Contact> contacts = structure.chainsWithAminoAcids()
                .map(chain -> chain.getFeature(PLIPInteractionContainer.class))
                .map(PLIPInteractionContainer::getInteractions)
                .flatMap(Collection::stream)
                // ignore covalently bound groups
                .filter(plipInteraction -> ContactsConstants.areNonCovalentGroups(plipInteraction.getPartner1(), plipInteraction.getPartner2()))
                .map(plipInteraction -> new Contact(((AminoAcid) plipInteraction.getPartner1()),
                        ((AminoAcid) plipInteraction.getPartner2()),
                        1.0,
                        plipInteraction.getClass().getSimpleName()))
                .collect(Collectors.toList());
        return new ContactMap(aminoAcids, contacts);
    }

    public static ContactMap createBetaCarbon8AngstromContactMap(Structure structure) {
        List<AminoAcid> aminoAcids = structure.aminoAcids().collect(Collectors.toList());
        List<Contact> contacts = ContactsConstants.determineNaiveInteractions(aminoAcids)
                .stream()
                .map(edge -> new Contact(edge.getLeft(), edge.getRight(), edge.getWeight(), "beta-carbon"))
                .collect(Collectors.toList());
        return new ContactMap(aminoAcids, contacts);
    }

    public List<AminoAcid> getAminoAcids() {
        return aminoAcids;
    }

    public List<Contact> getContacts() {
        return contacts;
    }

    public boolean hasContact(Pair<AminoAcid, AminoAcid> contact) {
        return hasContact(contact.getLeft(), contact.getRight());
    }

    public boolean hasContact(AminoAcid aminoAcid1, AminoAcid aminoAcid2) {
        return contacts.stream()
                .anyMatch(contact -> (contact.getLeft().equals(aminoAcid1) && contact.getRight().equals(aminoAcid2)) ||
                        (contact.getLeft().equals(aminoAcid2) && contact.getRight().equals(aminoAcid1)));
    }
}
