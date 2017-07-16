package studies.membrane.pdbtm.t02.fragments;

import de.bioforscher.jstructure.align.AlignmentPolicy;
import de.bioforscher.jstructure.align.StructureAlignmentQuery;
import de.bioforscher.jstructure.align.StructureAlignmentResult;
import de.bioforscher.jstructure.align.impl.SingleValueDecompositionAligner;
import de.bioforscher.jstructure.feature.interactions.PLIPInteraction;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.motif.SequenceMotif;
import de.bioforscher.jstructure.feature.motif.SequenceMotifDefinition;
import de.bioforscher.jstructure.model.identifier.IdentifierFactory;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Structure;
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

/**
 * Aim:
 * Extract all purely transmembrane fragments from the dataset.
 *
 * Results:
 * Directories of aligned sequence motifs with additional interaction data (represented by pseudo-atoms). Need to be
 * filtered (some observations are degenerated and do not really match the population). Furthermore, investigate where
 * certain interactions occur - ideally residue <code>i</code> should be interacting with residue <code>i + 4</code>,
 * what happens when the patterns moves to neighboring residues? What governs interactions on certain groups to be
 * absent?
 *
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
        Path tsvPath = outputPath.resolve("interactions-" + sequenceMotif.getMotifDefinition().name() + ".tsv");
        // create dirs if needed
        MembraneConstants.createDirectories(sequenceMotifDirectory);
        // write tsv file if needed
        if(!Files.exists(tsvPath)) {
            MembraneConstants.write(tsvPath, "id\ttype\toriginPos\ttargetPos" + System.lineSeparator());
        }

        // determine output filename
        AminoAcid firstResidue = sequenceMotif.getAminoAcids().get(0);
        AminoAcid lastResidue = sequenceMotif.getAminoAcids().get(sequenceMotif.getAminoAcids().size() - 1);
        String id = sequenceMotif.getChainId().getFullName();
        String filename = id + "-" + firstResidue.getIdentifier() + "-" + lastResidue.getIdentifier() + ".pdb";
        Path sequenceMotifOutputPath = sequenceMotifDirectory.resolve(filename);

        // if map does not contain reference already: make the current sequence motif reference
        Structure fragment = sequenceMotif.getAminoAcids().stream()
                .collect(StructureCollectors.toIsolatedStructure());
        fragment.setIdentifier(id);
        if(!referenceMotifs.containsKey(sequenceMotif.getMotifDefinition())) {
            referenceMotifs.put(sequenceMotif.getMotifDefinition(), fragment);
        }

        // align container
        GroupContainer reference = referenceMotifs.get(sequenceMotif.getMotifDefinition());
        logger.info("[{}] aligning {} to reference {}",
                sequenceMotif.getMotifDefinition().name(),
                filename.split("\\.")[0],
                reference.getIdentifier() + "-" +
                        reference.getAtoms().get(0).getParentGroup().getIdentifier() + "-" +
                        reference.getAtoms().get(reference.getAtoms().size() - 1).getParentGroup().getIdentifier());

        StructureAlignmentQuery query = StructureAlignmentQuery.of(reference, fragment)
                .matchingBehavior(AlignmentPolicy.MatchingBehavior.aminoAcidsComparableBackboneAtomNames)
                .manipulationBehavior(AlignmentPolicy.ManipulationBehavior.COPY);
        StructureAlignmentResult alignment = new SingleValueDecompositionAligner().align(query);

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
                .filter(plipInteraction -> plipInteraction.getPartner1().getResidueIdentifier().getResidueNumber() <
                        plipInteraction.getPartner2().getResidueIdentifier().getResidueNumber())
                .collect(Collectors.toList());

        // store interactions in synthetic group
        Group interactionGroup = new Group("INT", IdentifierFactory.createResidueIdentifier(999), true);
        interactions.stream()
                .map(PLIPInteraction::getAtomRepresentation)
                .forEach(interactionGroup::addAtom);
        // employ transformation on whole fragment including interactions
        GroupContainer alignedFragment = fragment.createDeepCopy();
        alignment.getTransformation().transform(alignedFragment);
        alignment.getTransformation().transform(interactionGroup);

        // write fragment file
        MembraneConstants.write(sequenceMotifOutputPath, alignedFragment.getPdbRepresentation() +
                interactionGroup.getPdbRepresentation());

        // append total count file
        try {
            String lines = interactions.stream()
                    .map(interaction -> id + "-" + firstResidue.getIdentifier() + "-" + lastResidue.getIdentifier() + "\t" +
                            interaction.getAtomRepresentation().getName() + "\t" +
                            sequenceMotif.getAminoAcids().indexOf(interaction.getPartner1()) + "\t" +
                                    sequenceMotif.getAminoAcids().indexOf(interaction.getPartner2()) +
                            System.lineSeparator())
                    .distinct()
                    .collect(Collectors.joining());

            // keep track of empty records
            if(lines.isEmpty()) {
                lines = filename.split("\\.")[0] + "\t-" + System.lineSeparator();
            }

            Files.write(tsvPath, lines.getBytes(), StandardOpenOption.APPEND);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }
}
