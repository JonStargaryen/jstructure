package de.bioforscher.jstructure.feature.uniprot.homologous;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry;

import java.util.List;

/**
 * Chain-specific mapping to homologous UniProt entries.
 * Created by bittrich on 7/10/17.
 */
public class UniProtHomologousEntryContainer extends FeatureContainerEntry {
    private final List<UniProtEntry> uniProtEntries;

    public UniProtHomologousEntryContainer(AbstractFeatureProvider featureProvider, List<UniProtEntry> uniProtEntries) {
        super(featureProvider);
        this.uniProtEntries = uniProtEntries;
    }

    public List<UniProtEntry> getUniProtEntries() {
        return uniProtEntries;
    }
}
