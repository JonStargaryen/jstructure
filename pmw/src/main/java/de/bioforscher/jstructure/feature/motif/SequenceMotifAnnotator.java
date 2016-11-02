package de.bioforscher.jstructure.feature.motif;

import de.bioforscher.jstructure.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.GroupContainer;
import de.bioforscher.jstructure.model.structure.Residue;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;

/**
 * Searches for sequence fragments matching any {@link SequenceMotifDefinition} and reports hits.
 * Created by S on 02.10.2016.
 */
public class SequenceMotifAnnotator implements FeatureProvider<GroupContainer> {
    private static final Logger LOGGER = Logger.getLogger(SequenceMotifAnnotator.class.getName());

    public enum FeatureNames {
        SEQUENCE_MOTIF
    }

    @Override
    @SuppressWarnings("unchecked")
    public void process(GroupContainer residueContainer) {
        // we need chain information as motifs cannot be part of 2 chains
        Map<String, List<Residue>> chains = residueContainer.residues().collect(Collectors.groupingBy(residue ->
                residue.getParentChain().getChainId()));
        for (String chainId : chains.keySet()) {
            List<Residue> residueInChain = chains.get(chainId);
            int chainLength = residueInChain.size();
            for (int resNum = 0; resNum < chainLength; resNum++) {
                Residue startResidue = residueInChain.get(resNum);
                // even though we express 1-letter-codes as Strings, excessive matching is probably much faster when done with chars
                char startAminoAcid = startResidue.getAminoAcid().getOneLetterCode().charAt(0);
                for (SequenceMotifDefinition candidate : SequenceMotifDefinition.values()) {
                    // get motif length
                    int motifLength = Integer.parseInt(candidate.name().substring(2));
                    char motifStart = candidate.name().charAt(0);
                    char motifEnd = candidate.name().charAt(1);

                    // chain not long enough to cover the proposed motif
                    if (resNum + motifLength >= chainLength) {
                        continue;
                    }

                    // start amino acid names does not match the motif
                    if (startAminoAcid != motifStart) {
                        continue;
                    }

                    Residue endResidue = residueInChain.get(resNum + motifLength);

                    // end amino acid does not match
                    if (endResidue.getAminoAcid().getOneLetterCode().charAt(0) != motifEnd) {
                        continue;
                    }

                    List<Residue> residueList = residueInChain.subList(resNum, resNum + motifLength + 1);
                    SequenceMotif sequenceMotif = new SequenceMotif(candidate, startResidue, endResidue);
                    residueList.forEach(residue -> {
                        List<SequenceMotif> value = residue.getFeature(List.class, FeatureNames.SEQUENCE_MOTIF);
                        // entry will be null at first - create list and assign reference
                        if(value == null) {
                            value = new ArrayList<>();
                            residue.setFeature(FeatureNames.SEQUENCE_MOTIF, value);
                        }
                        value.add(sequenceMotif);

                    });
                    LOGGER.log(Level.FINER, "found sequence motif: " + sequenceMotif);
                }
            }
        }
    }
}