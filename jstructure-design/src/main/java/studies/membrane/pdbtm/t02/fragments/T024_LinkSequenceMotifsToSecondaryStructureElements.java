package studies.membrane.pdbtm.t02.fragments;

import de.bioforscher.jstructure.feature.motif.SequenceMotif;
import de.bioforscher.jstructure.feature.motif.SequenceMotifContainer;
import de.bioforscher.jstructure.feature.sse.SecondaryStructureElement;
import de.bioforscher.jstructure.feature.sse.assp.ASSPSecondaryStructure;
import de.bioforscher.jstructure.feature.sse.assp.AssignmentOfSecondaryStructureInProteins;
import de.bioforscher.jstructure.feature.topology.Topology;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.ProteinParser;
import org.jsoup.Jsoup;
import studies.membrane.MembraneConstants;

import java.nio.file.Path;
import java.util.Collection;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * ASSP provides more fine-grained annotation of secondary structure elements especially with regards to helices. Use
 * this annotation to study sequence motifs. Are there preferences for certain sequence motifs to adapt non-standard
 * helix conformations?
 *
 * Results:
 * in {@link T024A_DistributionStatisticsOfAlphaThreePiHelices}
 *
 * Created by bittrich on 7/5/17.
 */
public class T024_LinkSequenceMotifsToSecondaryStructureElements {
    static final Path OUTPUT_PATH = MembraneConstants.PDBTM_FRAGMENTS_TM_BY_SEQUENCE_MOTIF_PATH.resolve("sse.tsv");

    public static void main(String[] args) {
        String output = getChainsAnnotatedByASSP()
                .map(Chain::getFeatureContainer)
                .map(featureContainer -> featureContainer.getFeature(SequenceMotifContainer.class))
                .map(SequenceMotifContainer::getSequenceMotifs)
                .flatMap(Collection::stream)
                // filter for tm
                .filter(sequenceMotif -> sequenceMotif.getAminoAcids()
                        .stream()
                        .allMatch(aminoAcid -> aminoAcid.getFeatureContainer().getFeature(Topology.class).isTransmembrane()))
                .map(T024_LinkSequenceMotifsToSecondaryStructureElements::handleSequenceMotif)
                .collect(Collectors.joining(System.lineSeparator(),
                        "id\tmotif\tsequence\t" +
//                                "dssp\t" +
                                "assp\tbend\tbendPos\tlength\talphaCount\tthreeCount\tpiCount\tcoilCount\t" +
                                "imperfect\tthree\tpi" + System.lineSeparator(),
                        ""));

        MembraneConstants.write(OUTPUT_PATH, output);
    }

    static Stream<Chain> getChainsAnnotatedByASSP() {
        return MembraneConstants.PdbtmAlphaNr.getIds()
                .map(T024_LinkSequenceMotifsToSecondaryStructureElements::mapToOptionalChain)
                .filter(Optional::isPresent)
                .map(Optional::get);
    }

    private static Optional<Chain> mapToOptionalChain(String id) {
        try {
            System.out.println("fetching " + id);
            String pdbId = id.split("_")[0];
            String chainId = id.split("_")[1];
            Protein protein = ProteinParser.source(MembraneConstants.PDBTM_PDB_PATH.resolve(pdbId + ".pdb"))
                    .minimalParsing(true)
                    .parse();
            Chain chain = protein.select()
                    .chainName(chainId)
                    .asChain();

            new AssignmentOfSecondaryStructureInProteins().process(protein);
//            MembraneConstants.SECONDARY_STRUCTURE_ANNOTATOR.process(protein);
            MembraneConstants.ORIENTATIONS_OF_PROTEINS_IN_MEMBRANES_ANNOTATOR.process(protein, Jsoup.parse(MembraneConstants.PDBTM_OPM_PATH.resolve(pdbId + ".xml").toFile(), "UTF-8"));
            MembraneConstants.SEQUENCE_MOTIF_ANNOTATOR.process(protein);

            return Optional.of(chain);
        } catch (Exception e) {
            // if computation fails or any file is missing
            return Optional.empty();
        }
    }

    private static String handleSequenceMotif(SequenceMotif sequenceMotif) {
//        List<DSSPSecondaryStructure> dsspSecondaryStructures = sequenceMotif.getAminoAcids()
//                .stream()
//                .map(AminoAcid::getFeatureContainer)
//                .map(featureContainer -> featureContainer.getFeature(DSSPSecondaryStructure.class))
//                .collect(Collectors.toList());
        List<ASSPSecondaryStructure> asspSecondaryStructures = sequenceMotif.getAminoAcids()
                .stream()
                .map(AminoAcid::getFeatureContainer)
                .map(featureContainer -> featureContainer.getFeature(ASSPSecondaryStructure.class))
                .collect(Collectors.toList());

        List<Integer> bends = asspSecondaryStructures.stream()
                .filter(asspSecondaryStructure -> asspSecondaryStructure.getBendAngle() > 60)
                .map(asspSecondaryStructures::indexOf)
                .mapToInt(i -> i + 1)
                .boxed()
                .collect(Collectors.toList());
        int length = asspSecondaryStructures.size();
        int alphaCount = (int) asspSecondaryStructures.stream().filter(asspSecondaryStructure -> asspSecondaryStructure.getSecondaryStructure() == SecondaryStructureElement.ALPHA_HELIX).count();
        int threeCount = (int) asspSecondaryStructures.stream().filter(asspSecondaryStructure -> asspSecondaryStructure.getSecondaryStructure() == SecondaryStructureElement.THREE_TEN_HELIX).count();
        int piCount = (int) asspSecondaryStructures.stream().filter(asspSecondaryStructure -> asspSecondaryStructure.getSecondaryStructure() == SecondaryStructureElement.PI_HELIX).count();
        int coilCount = (int) asspSecondaryStructures.stream().filter(asspSecondaryStructure -> asspSecondaryStructure.getSecondaryStructure() == SecondaryStructureElement.COIL).count();

        String output = sequenceMotif.getChainId() + "_" + sequenceMotif.getStartResidueNumber() + "-" + sequenceMotif.getEndResidueNumber() + "\t" +
                sequenceMotif.getMotifDefinition() + "\t" +
                sequenceMotif.getAminoAcids().stream().map(AminoAcid::getOneLetterCode).collect(Collectors.joining()) + "\t" +
//                dsspSecondaryStructures.stream().map(DSSPSecondaryStructure::getSecondaryStructure).map(SecondaryStructureElement::getOneLetterRepresentation).collect(Collectors.joining()) + "\t" +
                asspSecondaryStructures.stream().map(ASSPSecondaryStructure::getSecondaryStructure).map(SecondaryStructureElement::getOneLetterRepresentation).collect(Collectors.joining()) + "\t" +
                !(bends.isEmpty()) + "\t" +
                bends + "\t" +
                length + "\t" +
                alphaCount + "\t" +
                threeCount + "\t" +
                piCount + "\t" +
                coilCount + "\t" +
                (alphaCount != length) + "\t" +
                (threeCount > 0) + "\t" +
                (piCount > 0);

        System.out.println(output);

        return output;
    }
}
