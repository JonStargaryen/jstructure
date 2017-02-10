package de.bioforscher.jstructure.parser.plip;

import de.bioforscher.jstructure.model.structure.Group;

/**
 * Describes one observation of a PLIP interaction.
 * Created by bittrich on 2/9/17.
 */
public class PLIPInteraction {
    private final PLIPInteractionType plipInteractionType;
    private final Group partner1;
    private final Group partner2;

    PLIPInteraction(PLIPInteractionType plipInteractionType, Group group1, Group group2) {
        //TODO consider interacting atoms?
        this.plipInteractionType = plipInteractionType;
        this.partner1 = group1;
        this.partner2 = group2;
    }

    public PLIPInteractionType getPlipInteractionType() {
        return plipInteractionType;
    }

    public String getPlipInteractionShortName() {
        return plipInteractionType.getShortName();
    }

    public Group getPartner1() {
        return partner1;
    }

    public Group getPartner2() {
        return partner2;
    }

    @Override
    public String toString() {
        return getClass().getSimpleName() + " type='" + plipInteractionType + "' partner1='" + partner1.getIdentifier() + "' partner2='" + partner2.getIdentifier() + "'";
    }
}
