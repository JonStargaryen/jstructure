package studies.gmlvq.fingerprint;

import de.bioforscher.jstructure.model.structure.identifier.ChainIdentifier;
import de.bioforscher.jstructure.model.structure.identifier.ProteinIdentifier;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Composes data sets for the finger print data.
 * Created by bittrich on 6/23/17.
 */
public class FingerPrintDataSetComposer {
    private static final String DECOY = "decoy";
    private static final String FINGERPRINT = "fingerprint";

    /**
     *
     * @param motifPath the path to the motif directory - supposed to contain 2 subdirectories: decoy and fingerprint
     */
    public FingerPrintDataSetComposer(Path motifPath) throws IOException {
        Path decoyDirectory = motifPath.resolve(DECOY);
        Path fingerprintDirectory = motifPath.resolve(FINGERPRINT);

        Set<ChainIdentifier> decoyChains = getChainIdentifiers(decoyDirectory);
        Set<ChainIdentifier> fingerprintChains = getChainIdentifiers(fingerprintDirectory);


    }

    private Set<ChainIdentifier> getChainIdentifiers(Path directory) throws IOException {
        return Files.list(directory)
                .map(Path::toFile)
                .map(File::getName)
                .map(name -> name.split("-")[0])
                .map(id -> id.split("_"))
                .map(split -> ChainIdentifier.createFromChainId(ProteinIdentifier.createFromPdbId(split[0]), split[1]))
                .collect(Collectors.toSet());
    }
}
