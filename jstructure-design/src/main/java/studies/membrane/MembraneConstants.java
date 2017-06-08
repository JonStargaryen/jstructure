package studies.membrane;

import de.bioforscher.jstructure.feature.interactions.PLIPAnnotator;
import de.bioforscher.jstructure.feature.motif.SequenceMotif;
import de.bioforscher.jstructure.feature.motif.SequenceMotifAnnotator;
import de.bioforscher.jstructure.feature.motif.SequenceMotifContainer;
import de.bioforscher.jstructure.feature.sse.SecondaryStructureAnnotator;
import de.bioforscher.jstructure.feature.topology.MembraneContainer;
import de.bioforscher.jstructure.feature.topology.OrientationsOfProteinsInMembranesAnnotator;
import de.bioforscher.jstructure.feature.topology.Topology;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.jsoup.Jsoup;
import studies.StudyConstants;

import java.nio.file.Path;
import java.util.Collection;
import java.util.Optional;
import java.util.stream.Stream;

/**
 * Constants and functions shared by the membrane project.
 * Created by bittrich on 6/7/17.
 */
public class MembraneConstants extends StudyConstants {
    public static final Path PHD = GIT.resolve("phd_sb_repo");
    public static final Path PDBTM_ALPHA_NR_LIST = PHD.resolve("data/pdbtm_alpha_nr.list.txt");
    public static final Path PDBTM = PHD.resolve("data/pdbtm_alpha_nr/");
    public static final Path PDBTM_PDB_PATH = PDBTM.resolve("pdb");
    public static final Path PDBTM_CIF_PATH = PDBTM.resolve("cif");
    public static final Path PDBTM_OPM_PATH = PDBTM.resolve("opm");
    public static final Path PDBTM_PLIP_PATH = PDBTM.resolve("plip");
    public static final Path PDBTM_STATISTICS_PATH = PDBTM.resolve("statistics");
    private static final SecondaryStructureAnnotator SECONDARY_STRUCTURE_ANNOTATOR = new SecondaryStructureAnnotator();
    private static final PLIPAnnotator PLIP_ANNOTATOR = new PLIPAnnotator();
    private static final OrientationsOfProteinsInMembranesAnnotator ORIENTATIONS_OF_PROTEINS_IN_MEMBRANES_ANNOTATOR = new OrientationsOfProteinsInMembranesAnnotator();
    private static final SequenceMotifAnnotator SEQUENCE_MOTIF_ANNOTATOR = new SequenceMotifAnnotator();

    public static Stream<String> getIdsOfPdbtmAlphaNrList() {
        return lines(PDBTM_ALPHA_NR_LIST)
                // skip header and malformed lines
                .filter(line -> line.length() == 6 && !line.endsWith("_"));
    }

    public static Stream<Chain> getChainsOfPdbtmAlphaNrList() {
        return getIdsOfPdbtmAlphaNrList()
                .map(MembraneConstants::mapToOptionalChain)
                .filter(Optional::isPresent)
                .map(Optional::get);
    }

    private static Optional<Chain> mapToOptionalChain(String id) {
        try {
            System.out.println("fetching " + id);
            String pdbId = id.split("_")[0];
            String chainId = id.split("_")[1];
            Protein protein = ProteinParser.source(PDBTM_PDB_PATH.resolve(pdbId + ".pdb")).minimalParsing(true).parse();
            Chain chain = protein.select()
                    .chainName(chainId)
                    .asChain();

            // tied to residues
            SECONDARY_STRUCTURE_ANNOTATOR.process(protein);
            // assigns container directly to chain
            PLIP_ANNOTATOR.process(chain, Jsoup.parse(MembraneConstants.PDBTM_PLIP_PATH.resolve(chain.getChainId().getFullName() + ".xml").toFile(), "UTF-8"));
            // tied to the protein
            ORIENTATIONS_OF_PROTEINS_IN_MEMBRANES_ANNOTATOR.process(protein, Jsoup.parse(MembraneConstants.PDBTM_OPM_PATH.resolve(pdbId + ".xml").toFile(), "UTF-8"));
            // tied to protein
            SEQUENCE_MOTIF_ANNOTATOR.process(protein);

            return Optional.of(chain);
        } catch (Exception e) {
            // if computation fails or any file is missing
            return Optional.empty();
        }
    }

    public static Stream<AminoAcid> getAminoAcidsOfPdbtmAlphaNrList() {
        return getChainsOfPdbtmAlphaNrList()
                .flatMap(Chain::aminoAcids);
    }

    public static Stream<AminoAcid> getAminoAcidsOfPdbtmAlphaNrListTransmembrane() {
        return getAminoAcidsOfPdbtmAlphaNrList()
                .filter(aminoAcid -> aminoAcid.getParentChain()
                        .getParentProtein()
                        .getFeatureContainer()
                        .getFeature(MembraneContainer.class)
                        .isTransmembraneGroup(aminoAcid));
    }

    public static Stream<SequenceMotif> getSequenceMotifsOfPdbtmAlphaNrListTransmembrane() {
        return getChainsOfPdbtmAlphaNrList()
                .map(Chain::getFeatureContainer)
                .map(featureContainer -> featureContainer.getFeature(SequenceMotifContainer.class))
                .map(SequenceMotifContainer::getSequenceMotifs)
                .flatMap(Collection::stream)
                // filter for tm
                .filter(sequenceMotif -> sequenceMotif.getAminoAcids()
                        .stream()
                        .allMatch(aminoAcid -> aminoAcid.getFeatureContainer().getFeature(Topology.class).isTransmembrane()));
    }
}
