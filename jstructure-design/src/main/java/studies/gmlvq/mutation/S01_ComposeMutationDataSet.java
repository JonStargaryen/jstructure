package studies.gmlvq.mutation;

import de.bioforscher.jstructure.feature.uniprot.homologous.UniProtBlastQuery;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.ProteinParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import studies.StudyConstants;
import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry;
import uk.ac.ebi.kraken.interfaces.uniprot.features.FeatureType;

import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

/**
 * Annotate features
 * Created by bittrich on 7/10/17.
 */
public class S01_ComposeMutationDataSet {
    public static void main(String[] args) throws IOException {
        Files.lines(StudyConstants.MUTATION_DATASET_SCHAEFER2012)
                .filter(line -> !line.startsWith("ID,"))
                .limit(1)
                .forEach(S01_ComposeMutationDataSet::handleLine);
    }

    private static void handleLine(String line) {
        String[] split = line.split(",");

        String pdbId = split[0].substring(0, 4);
        String chainId = split[0].substring(4, 5);
        int position = Integer.valueOf(split[1]);
        AminoAcid.Family mutationTarget = AminoAcid.Family.resolveOneLetterCode(split[2]);
        double ddG = Double.valueOf(split[3]);

        Chain chain = ProteinParser.source(pdbId).parse()
                .select()
                .chainName(chainId)
                .asChain();

        List<UniProtEntry> hits = new UniProtBlastQuery().runUniProtBlastService(chain.getAminoAcidSequence());
        System.out.println(chain.getChainId() + " " + hits.size() + " similar sequences");
        hits.stream()
                .map(entry -> entry.getPrimaryUniProtAccession() + " " + entry.getFeatures(FeatureType.MUTAGEN).size() +
                        " mutations " + entry.getFeatures(FeatureType.VARIANT).size() + " variants")
                .forEach(System.out::println);
    }
}
