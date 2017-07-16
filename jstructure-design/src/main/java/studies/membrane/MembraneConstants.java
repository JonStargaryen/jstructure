package studies.membrane;

import de.bioforscher.jstructure.feature.interactions.PLIPAnnotator;
import de.bioforscher.jstructure.feature.interactions.PLIPInteraction;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.motif.SequenceMotif;
import de.bioforscher.jstructure.feature.motif.SequenceMotifAnnotator;
import de.bioforscher.jstructure.feature.motif.SequenceMotifContainer;
import de.bioforscher.jstructure.feature.sse.dssp.DictionaryOfProteinSecondaryStructure;
import de.bioforscher.jstructure.feature.topology.MembraneContainer;
import de.bioforscher.jstructure.feature.topology.OrientationsOfProteinsInMembranesAnnotator;
import de.bioforscher.jstructure.feature.topology.Topology;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
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
    public static final Path PDBTM_ALPHA_NR_LIST = PHD_REPO.resolve("data/pdbtm_alpha_nr.list.txt");
    public static final Path PDBTM = PHD_REPO.resolve("data/pdbtm_alpha_nr/");
    public static final Path PDBTM_PDB_PATH = PDBTM.resolve("pdb");
    public static final Path PDBTM_CIF_PATH = PDBTM.resolve("cif");
    public static final Path PDBTM_OPM_PATH = PDBTM.resolve("opm");
    public static final Path PDBTM_PLIP_PATH = PDBTM.resolve("plip");
    public static final Path PDBTM_STATISTICS_PATH = PDBTM.resolve("statistics");
    public static final Path PDBTM_FRAGMENTS_PATH = PDBTM.resolve("fragments");
    public static final Path PDBTM_FRAGMENTS_TM_PATH = PDBTM_FRAGMENTS_PATH.resolve("tm");
    public static final Path PDBTM_FRAGMENTS_TM_BY_SEQUENCE_MOTIF_PATH = PDBTM_FRAGMENTS_TM_PATH.resolve("by_sequenceMotif");
    public static final DictionaryOfProteinSecondaryStructure SECONDARY_STRUCTURE_ANNOTATOR = new DictionaryOfProteinSecondaryStructure();
    public static final PLIPAnnotator PLIP_ANNOTATOR = new PLIPAnnotator();
    public static final OrientationsOfProteinsInMembranesAnnotator ORIENTATIONS_OF_PROTEINS_IN_MEMBRANES_ANNOTATOR = new OrientationsOfProteinsInMembranesAnnotator();
    public static final SequenceMotifAnnotator SEQUENCE_MOTIF_ANNOTATOR = new SequenceMotifAnnotator();

    private static Optional<Chain> mapToOptionalChain(String id) {
        try {
            System.out.println("fetching " + id);
            String pdbId = id.split("_")[0];
            String chainId = id.split("_")[1];
            Structure protein = StructureParser.source(PDBTM_PDB_PATH.resolve(pdbId + ".pdb"))
                    .minimalParsing(true)
                    .parse();
            Chain chain = protein.select()
                    .chainName(chainId)
                    .asChain();

            // tied to residues
            SECONDARY_STRUCTURE_ANNOTATOR.process(protein);
            // assigns container directly to chain
            PLIP_ANNOTATOR.process(chain, Jsoup.parse(MembraneConstants.PDBTM_PLIP_PATH.resolve(chain.getChainIdentifier().getFullName() + ".xml").toFile(), "UTF-8"));
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

    public static class PdbtmAlphaNr {
        public static Stream<String> getIds() {
            return lines(PDBTM_ALPHA_NR_LIST)
                    // skip header and malformed lines
                    .filter(line -> line.length() == 6 && !line.endsWith("_"));
        }

        public static Stream<Chain> getChains() {
            return getIds()
                    .map(MembraneConstants::mapToOptionalChain)
                    .filter(Optional::isPresent)
                    .map(Optional::get);
        }

        public static Stream<AminoAcid> getAminoAcids() {
            return getChains()
                    .flatMap(Chain::aminoAcids);
        }

        public static Stream<AminoAcid> getAminoAcidsTransmembrane() {
            return getAminoAcids()
                    .filter(aminoAcid -> aminoAcid.getParentChain()
                            .getParentStructure()
                            .getFeatureContainer()
                            .getFeature(MembraneContainer.class)
                            .isTransmembraneGroup(aminoAcid));
        }

        public static Stream<SequenceMotif> getSequenceMotifs() {
            return getChains()
                    .map(Chain::getFeatureContainer)
                    .map(featureContainer -> featureContainer.getFeature(SequenceMotifContainer.class))
                    .map(SequenceMotifContainer::getSequenceMotifs)
                    .flatMap(Collection::stream);
        }

        public static Stream<SequenceMotif> getSequenceMotifsTransmembrane() {
            return getSequenceMotifs()
                    // filter for tm
                    .filter(sequenceMotif -> sequenceMotif.getAminoAcids()
                            .stream()
                            .allMatch(aminoAcid -> aminoAcid.getFeatureContainer().getFeature(Topology.class).isTransmembrane()));
        }

        public static Stream<PLIPInteraction> getInteractions() {
            return getChains()
                    .map(Chain::getFeatureContainer)
                    .map(featureContainer -> featureContainer.getFeature(PLIPInteractionContainer.class))
                    .map(PLIPInteractionContainer::getInteractions)
                    .flatMap(Collection::stream);
        }

        public static Stream<PLIPInteraction> getInteractionsTransmembrane() {
            return getInteractions()
                    .filter(plipInteraction -> plipInteraction.getPartner1().isAminoAcid() && plipInteraction.getPartner2().isAminoAcid())
                    .filter(plipInteraction -> plipInteraction.getPartner1().getFeatureContainer().getFeature(Topology.class).isTransmembrane() &&
                        plipInteraction.getPartner2().getFeatureContainer().getFeature(Topology.class).isTransmembrane());
        }

        public static Stream<PLIPInteraction> getInteractionsBetweenTransmembraneHelices() {
            return getInteractionsTransmembrane()
                    // filter for interactions of different parts of transmembrane helices
                    .filter(plipInteraction -> {
                        MembraneContainer membraneContainer = plipInteraction.getPartner1().getParentChain().getParentStructure().getFeatureContainer().getFeature(MembraneContainer.class);
                        return membraneContainer.getEmbeddingTransmembraneSegment(plipInteraction.getPartner1()) != membraneContainer.getEmbeddingTransmembraneSegment(plipInteraction.getPartner2());
                    });
        }
    }
}
