package de.bioforscher.jstructure.graph.contact.definition;

public class ContactDefinitionFactory {
    private ContactDefinitionFactory() {
        // deny instantiation
    }

    public static ContactDefinition createAlphaCarbonContactDefinition(double distance) {
        return new AlphaCarbonContactDefinition(distance);
    }

    public static ContactDefinition createBetaCarbonContactDefinition(double distance) {
        return new BetaCarbonContactDefinition(distance);
    }

    public static ContactDefinition createHydrogenBondContactDefinition() {
        return new HydrogenBondContactDefinition();
    }

    public static ContactDefinition createHydrophobicInteractionContactDefinition() {
        return new HydrophobicInteractionContactDefinition();
    }

    public static ContactDefinition createPlipContactDefinition() {
        return new PlipContactDefinition();
    }
}
