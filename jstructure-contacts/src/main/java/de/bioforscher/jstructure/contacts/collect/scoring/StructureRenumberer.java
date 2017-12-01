package de.bioforscher.jstructure.contacts.collect.scoring;

import de.bioforscher.jstructure.model.identifier.IdentifierFactory;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.List;
import java.util.stream.Collectors;

public class StructureRenumberer {
    private static final Logger logger = LoggerFactory.getLogger(StructureRenumberer.class);

    public static void renumberStructure(Structure reference, Structure model) {
        List<AminoAcid> modelAminoAcids = model.aminoAcids().collect(Collectors.toList());
        String referenceSequence = reference.getAminoAcidSequence();
        String modelSequence = model.getAminoAcidSequence();
        int referenceLength = referenceSequence.length();
        int modelLength = modelSequence.length();
        logger.info("renumbering model with length {} wrt to reference with length {}",
                referenceLength,
                modelLength);

        try {
            SequencePair<ProteinSequence, AminoAcidCompound> pair = Alignments.getPairwiseAlignment(new ProteinSequence(referenceSequence),
                    new ProteinSequence(modelSequence),
                    Alignments.PairwiseSequenceAlignerType.GLOBAL,
                    new SimpleGapPenalty(),
                    SubstitutionMatrixHelper.getBlosum62());

            logger.info("alignment:{}{}",
                    System.lineSeparator(),
                    pair.toString());

            String alignedModelSequence = pair.getQuery().getSequenceAsString();
            int sequencePosition = 0;
            for(AminoAcid aminoAcid : modelAminoAcids) {
                String oldResidueNumber = aminoAcid.toString();
                while(sequencePosition < alignedModelSequence.length() && alignedModelSequence.charAt(sequencePosition) == '-') {
                    sequencePosition++;
                }
                ResidueIdentifier newResidueIdentifier = IdentifierFactory.createResidueIdentifier(sequencePosition + 1);
                aminoAcid.setResidueIdentifier(newResidueIdentifier);
                logger.debug("renumbered {} to {}",
                        oldResidueNumber,
                        aminoAcid.toString());
                sequencePosition++;
            }
        } catch (CompoundNotFoundException e) {
            throw new IllegalArgumentException(e);
        }
    }
}
