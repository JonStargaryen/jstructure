package studies.membrane.pdbtm.t02.fragments;

import de.bioforscher.jstructure.feature.motif.SequenceMotif;
import de.bioforscher.jstructure.feature.motif.SequenceMotifContainer;
import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.feature.sse.SecondaryStructureElement;
import de.bioforscher.jstructure.feature.sse.assp.AssignmentOfSecondaryStructureInProteinsWrapper;
import de.bioforscher.jstructure.feature.topology.Topology;
import de.bioforscher.jstructure.model.Combinatorics;
import de.bioforscher.jstructure.model.Fragment;
import de.bioforscher.jstructure.model.feature.AbstractFeatureable;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.ProteinParser;
import org.jsoup.Jsoup;
import studies.membrane.MembraneConstants;

import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Look at fragments
 * Created by bittrich on 7/5/17.
 */
public class T025_HelixDistribution {
    public static void main(String[] args) {
        String output = Stream.concat(getTransmembraneThreeTenHelices(), getTransmembranePiHelices())
                .map(T025_HelixDistribution::composeLine)
                .peek(System.out::println)
                .collect(Collectors.joining(System.lineSeparator(),
                        "id\tsse\tsequence\tmotifs1\tmotifs2\tmotifs3\tmotifs4\tmotifs5" + System.lineSeparator(),
                        ""));

        System.out.println(output);
    }

    private static String composeLine(Fragment<AminoAcid> aminoAcidFragment) {
        String id = aminoAcidFragment.getElement(0).getParentChain().getChainId() + "-" + aminoAcidFragment.getElement(0).getResidueIdentifier().getResidueNumber() + "-" + aminoAcidFragment.getElement(aminoAcidFragment.getElements().size() - 1).getResidueIdentifier().getResidueNumber();
        String secondaryStructure = aminoAcidFragment.getElement(0).getFeatureContainer().getFeature(GenericSecondaryStructure.class).getSecondaryStructure().name();
        String sequence = aminoAcidFragment.getElements()
                .stream()
                .map(AminoAcid::getOneLetterCode)
                .collect(Collectors.joining());
        SequenceMotifContainer sequenceMotifContainer = aminoAcidFragment.getElement(0).getParentChain().getFeatureContainer().getFeature(SequenceMotifContainer.class);
        List<List<String>> sequenceMotifs = aminoAcidFragment.getElements()
                .stream()
                .map(sequenceMotifContainer::getEmbeddingSequenceMotifsFor)
                .map(container -> container.getSequenceMotifs()
                        .stream()
                        .map(SequenceMotif::getMotifDefinition)
                        .map(SequenceMotifDefinition::name)
                        .collect(Collectors.toList()))
                .collect(Collectors.toList());
        return id + "\t" + secondaryStructure + "\t" +
                sequence + "\t" +
                sequenceMotifs.stream()
                        .map(Object::toString)
                        .collect(Collectors.joining("\t"));
    }

    static Stream<Fragment<AminoAcid>> getTransmembranePiHelices() {
        return getChainsAnnotatedByASSP()
                // create all 5mers of this chain
                .flatMap(chain -> {
                    try {
                        return Combinatorics.fragmentsOf(chain.aminoAcids()
                                .filter(aminoAcid -> aminoAcid.getFeatureContainer().getFeature(Topology.class).isTransmembrane())
                                .collect(Collectors.toList()), 5);
                    } catch (IllegalArgumentException e) {
                        return Stream.empty();
                    }
                })
                // ensure consecutive fragments
                .filter(fragment -> IntStream.range(1, fragment.getElements().size())
                        .allMatch(i -> fragment.getElement(i).getResidueIdentifier().getResidueNumber() == fragment.getElement(i - 1).getResidueIdentifier().getResidueNumber() + 1))
                // filter for fragments exclusively
                .filter(fragment -> fragment.getElements().stream()
                        .map(AbstractFeatureable::getFeatureContainer)
                        .map(featureContainer -> featureContainer.getFeature(GenericSecondaryStructure.class))
                        .map(GenericSecondaryStructure::getSecondaryStructure)
                        .allMatch(secondaryStructureElement -> secondaryStructureElement == SecondaryStructureElement.PI_HELIX));
    }

    private static Stream<Chain> getChainsAnnotatedByASSP() {
        return MembraneConstants.PdbtmAlphaNr.getIds()
                    .map(T025_HelixDistribution::mapToOptionalChain)
                    .filter(Optional::isPresent)
                    .map(Optional::get);
    }

    private static Optional<Chain> mapToOptionalChain(String id) {
        try {
            System.out.println("processing " + id);
            String pdbId = id.split("_")[0];
            String chainId = id.split("_")[1];
            Protein protein = ProteinParser.source(MembraneConstants.PDBTM_PDB_PATH.resolve(pdbId + ".pdb"))
                    .minimalParsing(true)
                    .parse();
            Chain chain = protein.select()
                    .chainName(chainId)
                    .asChain();

            new AssignmentOfSecondaryStructureInProteinsWrapper().process(protein);
            MembraneConstants.ORIENTATIONS_OF_PROTEINS_IN_MEMBRANES_ANNOTATOR.process(protein, Jsoup.parse(MembraneConstants.PDBTM_OPM_PATH.resolve(pdbId + ".xml").toFile(), "UTF-8"));
            MembraneConstants.SEQUENCE_MOTIF_ANNOTATOR.process(protein);

            return Optional.of(chain);
        } catch (Exception e) {
            // if computation fails or any file is missing
            return Optional.empty();
        }
    }

    static Stream<Fragment<AminoAcid>> getTransmembraneThreeTenHelices() {
        return getChainsAnnotatedByASSP()
                // create all 5mers of this chain
                .flatMap(chain -> {
                    try {
                        return Combinatorics.fragmentsOf(chain.aminoAcids()
                                .filter(aminoAcid -> aminoAcid.getFeatureContainer().getFeature(Topology.class).isTransmembrane())
                                .collect(Collectors.toList()), 3);
                    } catch (IllegalArgumentException e) {
                        return Stream.empty();
                    }
                })
                // ensure consecutive fragments
                .filter(fragment -> IntStream.range(1, fragment.getElements().size())
                        .allMatch(i -> fragment.getElement(i).getResidueIdentifier().getResidueNumber() == fragment.getElement(i - 1).getResidueIdentifier().getResidueNumber() + 1))
                // filter for fragments exclusively
                .filter(fragment -> fragment.getElements().stream()
                        .map(AbstractFeatureable::getFeatureContainer)
                        .map(featureContainer -> featureContainer.getFeature(GenericSecondaryStructure.class))
                        .map(GenericSecondaryStructure::getSecondaryStructure)
                        .allMatch(secondaryStructureElement -> secondaryStructureElement == SecondaryStructureElement.THREE_TEN_HELIX));
    }
}
