package de.bioforscher.jstructure.feature.uniprot.homologous;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.structure.identifier.ChainIdentifier;
import de.bioforscher.jstructure.model.structure.identifier.IdentifierFactory;
import uk.ac.ebi.kraken.interfaces.uniprot.DatabaseType;
import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry;

import java.util.Collection;
import java.util.List;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * Chain-specific mapping to homologous UniProt entries.
 * Created by bittrich on 7/10/17.
 */
public class UniProtHomologousEntryContainer extends FeatureContainerEntry {
    private final List<UniProtEntry> uniProtEntries;
    private final List<ChainIdentifier> homologousChains;
    private static final Pattern CHAIN_PATTERN = Pattern.compile("/");

    public UniProtHomologousEntryContainer(AbstractFeatureProvider featureProvider, List<UniProtEntry> uniProtEntries) {
        super(featureProvider);
        this.uniProtEntries = uniProtEntries;
        this.homologousChains = uniProtEntries.stream()
                .map(uniProtEntry -> uniProtEntry.getDatabaseCrossReferences(DatabaseType.PDB))
                .flatMap(Collection::stream)
                .flatMap(databaseCrossReference -> {
                    //TODO range is also available in the format fourth -> "A=321-419"
                    String[] split = databaseCrossReference.getFourth().getValue().split("=");
                    // some entries reference multiple chains separated by '/'
                    return CHAIN_PATTERN.splitAsStream(split[0])
                            .map(chainId -> IdentifierFactory.createChainIdentifier(databaseCrossReference.getPrimaryId().getValue(), chainId));
                })
                .collect(Collectors.toList());
    }

    public List<UniProtEntry> getUniProtEntries() {
        return uniProtEntries;
    }

    public List<ChainIdentifier> getHomologousChains() {
        return homologousChains;
    }
}
