package de.bioforscher.jstructure.model.structure.prototype;

import de.bioforscher.jstructure.mathematics.LinearAlgebra;
import de.bioforscher.jstructure.model.feature.AbstractFeatureable;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.selection.Selection;

import java.util.List;

/**
 * Shared implementation of all {@link Group} implementations.
 * Created by bittrich on 5/24/17.
 */
public abstract class AbstractGroup extends AbstractFeatureable implements Group {
    private final String threeLetterCode;

    AbstractGroup(String threeLetterCode) {
        this.threeLetterCode = threeLetterCode;
    }

    @Override
    public String getThreeLetterCode() {
        return threeLetterCode;
    }

    static <G extends AbstractGroup> G createPrototypeInstance(String id) {

    }
}
