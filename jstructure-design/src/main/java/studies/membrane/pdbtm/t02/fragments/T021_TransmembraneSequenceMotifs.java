package studies.membrane.pdbtm.t02.fragments;

import de.bioforscher.jstructure.alignment.Alignment;
import de.bioforscher.jstructure.alignment.AlignmentPolicy;
import de.bioforscher.jstructure.alignment.StructureAligner;
import de.bioforscher.jstructure.feature.interactions.PLIPInteraction;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.motif.SequenceMotif;
import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.ResidueNumber;
import de.bioforscher.jstructure.model.structure.StructureCollectors;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import studies.membrane.MembraneConstants;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Extract all purely transmembrane fragments from the dataset.
 * Created by bittrich on 6/16/17.
 */
public class T021_TransmembraneSequenceMotifs {
    private static final Logger logger = LoggerFactory.getLogger(T021_TransmembraneSequenceMotifs.class);
    private static final Path outputPath = MembraneConstants.PDBTM_FRAGMENTS_TM_BY_SEQUENCE_MOTIF_PATH;
    private static final Map<SequenceMotifDefinition, GroupContainer> referenceMotifs = new HashMap<>();

    public static void main(String[] args) {
        purgeDirectory(outputPath.toFile());
        MembraneConstants.PdbtmAlphaNr.getSequenceMotifsTransmembrane()
                .forEach(T021_TransmembraneSequenceMotifs::handleSequenceMotif);
    }

    private static void purgeDirectory(File directory) {
        try {
            for (File file : directory.listFiles()) {
                if (file.isDirectory()) {
                    purgeDirectory(file);
                }
                file.delete();
            }
        } catch (NullPointerException e) {
            // nothing to purge, when directory does not exist yet
        }
    }

    private static void handleSequenceMotif(SequenceMotif sequenceMotif) {
        Path sequenceMotifDirectory = outputPath.resolve(sequenceMotif.getMotifDefinition().name());
        Path tsvPath = sequenceMotifDirectory.resolve("interactions.tsv");
        // create dirs if needed
        MembraneConstants.createDirectories(sequenceMotifDirectory);
        // write tsv file if needed
        if(!Files.exists(tsvPath)) {
            MembraneConstants.write(tsvPath, "id\ttype\toriginPos\ttargetPos" + System.lineSeparator());
        }

        // if map does not contain reference already: make the current sequence motif reference
        Group interactionGroup = new Group("INT", new ResidueNumber(999), true);
//        GroupContainer fragment = sequenceMotif.getAminoAcids().stream()
//                .collect(StructureCollectors.toGroupContainer());
        GroupContainer fragment = Stream.concat(sequenceMotif.getAminoAcids().stream(), Stream.of(interactionGroup))
                .collect(StructureCollectors.toGroupContainer());
        if(!referenceMotifs.containsKey(sequenceMotif.getMotifDefinition())) {
            referenceMotifs.put(sequenceMotif.getMotifDefinition(), fragment);
        }

        // determine output filename
        AminoAcid firstResidue = sequenceMotif.getAminoAcids().get(0);
        AminoAcid lastResidue = sequenceMotif.getAminoAcids().get(sequenceMotif.getAminoAcids().size() - 1);
        String id = firstResidue.getParentChain().getChainId().getFullName();
        String filename = id + "-" + firstResidue.getIdentifier() + "-" + lastResidue.getIdentifier() + ".pdb";
        Path sequenceMotifOutputPath = sequenceMotifDirectory.resolve(filename);

        List<Atom> sequenceMotifAtoms = fragment.atoms()
                .collect(Collectors.toList());
        List<PLIPInteraction> interactions = fragment.groups()
                .filter(AminoAcid.class::isInstance)
                .map(Group::getFeatureContainer)
                .map(featureContainer -> featureContainer.getFeature(PLIPInteractionContainer.class))
                .map(PLIPInteractionContainer::getInteractions)
                .flatMap(Collection::stream)
                // filter for pure-backbone interactions
                .filter(PLIPInteraction::isBackboneInteraction)
                // filter for interactions fully embedded into this sequence motif
                .filter(plipInteraction -> plipInteraction.allAtoms().allMatch(sequenceMotifAtoms::contains))
                // make list distinct
                .filter(plipInteraction -> plipInteraction.getPartner1().getResidueNumber().getResidueNumber() <
                        plipInteraction.getPartner2().getResidueNumber().getResidueNumber())
                .collect(Collectors.toList());
        interactions.stream()
                .map(PLIPInteraction::getAtomRepresentation)
                .forEach(interactionGroup::addAtom);

        // align container and all interactions to reference
        GroupContainer reference = referenceMotifs.get(sequenceMotif.getMotifDefinition());
        logger.info("[{}] aligning {} to reference {}",
                sequenceMotif.getMotifDefinition().name(),
                filename.split("\\.")[0],
                reference.getAtoms().get(0).getParentGroup().getParentChain().getChainId().getFullName() + "-" +
                        reference.getAtoms().get(0).getParentGroup().getIdentifier() + "-" +
                        reference.getAtoms().get(reference.getAtoms().size() - 1).getParentGroup().getIdentifier());

        Alignment alignment = StructureAligner.builder(reference, fragment)
                .matchingBehavior(AlignmentPolicy.MatchingBehavior.COMPARABLE_ATOM_NAMES)
                .manipulationBehavior(AlignmentPolicy.ManipulationBehavior.COPY)
                .align();
//        alignment.transform(interactionGroup);
        // write fragment file
        MembraneConstants.write(sequenceMotifOutputPath, alignment.getAlignedQuery().getPdbRepresentation() +
                interactionGroup.getPdbRepresentation());

        // append total count file
        try {
            String lines = interactions.stream()
                    .map(interaction -> id + "-" + firstResidue.getIdentifier() + "-" + lastResidue.getIdentifier() + "\t" +
                            sequenceMotif.getAminoAcids().indexOf(interaction.getPartner1()) + "\t" +
                                    sequenceMotif.getAminoAcids().indexOf(interaction.getPartner2()) +
                            System.lineSeparator())
                    .distinct()
                    .collect(Collectors.joining());
            Files.write(tsvPath, lines.getBytes(), StandardOpenOption.APPEND);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }
}
