package de.bioforscher.jstructure.efr.model;

import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;

import java.util.ArrayList;
import java.util.List;

public class Start2FoldResidueAnnotation extends FeatureContainerEntry {
    private final List<ProtectionLevel> protectionLevelEntries;

    public Start2FoldResidueAnnotation() {
        super(null);
        this.protectionLevelEntries = new ArrayList<>();
    }

    public List<ProtectionLevel> getProtectionLevelEntries() {
        return protectionLevelEntries;
    }

    public void addProtectionLevelEntry(ProtectionLevel protectionLevel) {
        protectionLevelEntries.add(protectionLevel);
    }

    public boolean isEarly() {
        return filter(ProtectionLevel.EARLY);
    }

    public boolean isIntermediate() {
        return filter(ProtectionLevel.INTERMEDIATE);
    }

    public boolean isLate() {
        return filter(ProtectionLevel.LATE);
    }

    public boolean isWeak() {
        return filter(ProtectionLevel.WEAK);
    }

    public boolean isMedium() {
        return filter(ProtectionLevel.MEDIUM);
    }

    public boolean isStrong() {
        return filter(ProtectionLevel.STRONG);
    }

    private boolean filter(ProtectionLevel protectionLevel) {
        return protectionLevelEntries.stream()
                .anyMatch(pl -> pl == protectionLevel);
    }
}