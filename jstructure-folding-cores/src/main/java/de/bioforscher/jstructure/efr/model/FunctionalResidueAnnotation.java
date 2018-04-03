package de.bioforscher.jstructure.efr.model;

import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;

import java.util.ArrayList;
import java.util.List;

public class FunctionalResidueAnnotation extends FeatureContainerEntry {
    private List<String> functionalAnnotations;

    public FunctionalResidueAnnotation() {
        super(null);
        this.functionalAnnotations = new ArrayList<>();
    }

    public List<String> getFunctionalAnnotations() {
        return functionalAnnotations;
    }

    public void addFunctionalAnnotation(String annotation) {
        functionalAnnotations.add(annotation);
    }

    public boolean isNonFunctional() {
        return functionalAnnotations.isEmpty();
    }

    public boolean isFunctional() {
        return !isNonFunctional();
    }
}