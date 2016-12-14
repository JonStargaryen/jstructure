package de.bioforscher.jstructure.feature.motif;

import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.AminoAcid;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Searches for sequence fragments matching any {@link SequenceMotifDefinition} and reports hits.
 * Created by S on 02.10.2016.
 */
public class SequenceMotifAnnotator implements FeatureProvider {
    private static final Logger logger = LoggerFactory.getLogger(SequenceMotifAnnotator.class);

    public enum FeatureNames {
        SEQUENCE_MOTIF
    }

    @Override
    @SuppressWarnings("unchecked")
    public void process(Protein protein) {
        List<SequenceMotif> globalListOfSequenceMotifs = new ArrayList<>();

        // we need getChain information as motifs cannot be part of 2 chains
        Map<String, List<Group>> chains = Selection.on(protein)
                .aminoAcids()
                .asFilteredGroups()
                .collect(Collectors.groupingBy(residue -> residue.getParentChain().getChainId()));
        for (String chainId : chains.keySet()) {
            List<Group> residueInChain = chains.get(chainId);
            int chainLength = residueInChain.size();
            for (int resNum = 0; resNum < chainLength; resNum++) {
                Group startResidue = residueInChain.get(resNum);
                // even though we express 1-letter-codes as Strings, excessive matching is probably much faster when done with chars
                char startAminoAcid = AminoAcid.valueOfIgnoreCase(startResidue.getPdbName()).getOneLetterCode().charAt(0);
                for (SequenceMotifDefinition candidate : SequenceMotifDefinition.values()) {
                    // findAny motif length
                    int motifLength = Integer.parseInt(candidate.name().substring(2));
                    char motifStart = candidate.name().charAt(0);
                    char motifEnd = candidate.name().charAt(1);

                    // getChain not long enough to cover the proposed motif
                    if (resNum + motifLength >= chainLength) {
                        continue;
                    }

                    // start amino acid names does not match the motif
                    if (startAminoAcid != motifStart) {
                        continue;
                    }

                    Group endResidue = residueInChain.get(resNum + motifLength);

                    // end amino acid does not match
                    if (AminoAcid.valueOfIgnoreCase(endResidue.getPdbName()).getOneLetterCode().charAt(0) != motifEnd) {
                        continue;
                    }

                    // motif lacks internal residues or has a invalid ordering
                    if(startResidue.getResidueNumber() + motifLength != endResidue.getResidueNumber()) {
                        continue;
                    }

                    List<Group> residueList = residueInChain.subList(resNum, resNum + motifLength + 1);
                    SequenceMotif sequenceMotif = new SequenceMotif(candidate, startResidue, endResidue);
                    residueList.forEach(residue -> {
                        List<SequenceMotif> value = residue.getFeature(List.class, FeatureNames.SEQUENCE_MOTIF);
                        // entry will be null at first - create list and assign reference
                        if(value == null) {
                            value = new ArrayList<>();
                            residue.setFeature(FeatureNames.SEQUENCE_MOTIF, value);
                        }
                        value.add(sequenceMotif);
                        if(!globalListOfSequenceMotifs.contains(sequenceMotif)) {
                            globalListOfSequenceMotifs.add(sequenceMotif);
                        }
                    });
                    logger.debug("found sequence motif: {}", sequenceMotif);
                }
            }
        }

        // set global reference to all entries
        protein.setFeature(FeatureNames.SEQUENCE_MOTIF, globalListOfSequenceMotifs);
    }
}