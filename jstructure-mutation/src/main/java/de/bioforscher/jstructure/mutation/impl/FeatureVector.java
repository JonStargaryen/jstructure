package de.bioforscher.jstructure.mutation.impl;

import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceArea;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfile;
import de.bioforscher.jstructure.feature.loopfraction.LoopFraction;
import de.bioforscher.jstructure.model.feature.FeatureContainer;
import de.bioforscher.jstructure.model.structure.Group;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.List;
import java.util.Locale;

/**
 * Describes properties on either a single amino acid or those of its environment.
 * Created by bittrich on 7/12/17.
 */
class FeatureVector {
    private static final DecimalFormat DECIMAL_FORMAT = new DecimalFormat("#.####", new DecimalFormatSymbols(Locale.US));
    private double rasa;
    private double energy;
    private double loopFraction;
    private double ligandContacts;
    //TODO more features

    FeatureVector(Group group) {
        FeatureContainer featureContainer = group.getFeatureContainer();
        this.rasa = featureContainer.getFeature(AccessibleSurfaceArea.class).getRelativeAccessibleSurfaceArea();
        this.energy = featureContainer.getFeature(EnergyProfile.class).getSolvationEnergy();
        this.loopFraction = featureContainer.getFeature(LoopFraction.class).getLoopFraction();
        this.ligandContacts = LigandContactScreener.determineNumberOfLigandContacts(group);
        //TODO more features
    }

    FeatureVector(List<Group> groups) {
        groups.stream()
                .map(FeatureVector::new)
                .forEach(this::consume);
        this.rasa /= groups.size();
        this.energy /= groups.size();
        this.loopFraction /= groups.size();
        this.ligandContacts /= groups.size();
        //TODO more features
    }

    /**
     * Constructs the delta vector of 2 feature vectors.
     * @param original
     * @param mutated
     */
    FeatureVector(FeatureVector original, FeatureVector mutated) {
        this.rasa = original.rasa - mutated.rasa;
        this.energy = original.energy - mutated.energy;
        this.loopFraction = original.loopFraction - mutated.loopFraction;
        this.ligandContacts = original.ligandContacts - mutated.ligandContacts;
        //TODO more features
    }

    private void consume(FeatureVector featureVector) {
        this.rasa += featureVector.rasa;
        this.energy += featureVector.energy;
        this.loopFraction += featureVector.loopFraction;
        this.ligandContacts += featureVector.ligandContacts;
        //TODO more features
    }

    public String toArffString() {
        return DECIMAL_FORMAT.format(rasa) + "," + DECIMAL_FORMAT.format(energy) + "," +
                DECIMAL_FORMAT.format(loopFraction) + "," + DECIMAL_FORMAT.format(ligandContacts);
    }
}
