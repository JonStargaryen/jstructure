package de.bioforscher.jstructure.feature.motif;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.family.GroupInformation;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Searches for sequence fragments matching any {@link SequenceMotifDefinition} and reports hits.
 * Created by S on 02.10.2016.
 */
@FeatureProvider(provides = SequenceMotifContainer.class)
public class SequenceMotifAnnotator extends AbstractFeatureProvider {
    private static final Logger logger = LoggerFactory.getLogger(SequenceMotifAnnotator.class);

    @Override
    protected void processInternally(Protein protein) {
        SequenceMotifContainer globalList = new SequenceMotifContainer(this);

        protein.aminoAcidChains()
                .forEach(chain -> {
            List<Group> aminoAcids = chain.aminoAcids().collect(Collectors.toList());

            // init with empty containers
            aminoAcids.forEach(group -> group.getFeatureContainer().addFeature(new SequenceMotifContainer(this)));

            int chainLength = aminoAcids.size();
            for (int resNum = 0; resNum < chainLength; resNum++) {
                Group startResidue = aminoAcids.get(resNum);
                // even though we express 1-letter-codes as Strings, excessive matching is probably much faster when done with chars
                char startAminoAcid = startResidue.getGroupInformation().getOneLetterCode().charAt(0);
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

                    Group endResidue = aminoAcids.get(resNum + motifLength);

                    // end amino acid does not match
                    if (endResidue.getGroupInformation().getOneLetterCode().charAt(0) != motifEnd) {
                        continue;
                    }

                    // motif lacks internal residues or has a invalid ordering
                    if(startResidue.getResidueNumber() + motifLength != endResidue.getResidueNumber()) {
                        continue;
                    }

                    List<Group> residueList = aminoAcids.subList(resNum, resNum + motifLength + 1);
                    SequenceMotif sequenceMotif = new SequenceMotif(candidate,
                            chain.getChainId(),
                            startResidue.getResidueNumber(),
                            endResidue.getResidueNumber(),
                            residueList.stream()
                                .map(Group::getGroupInformation)
                                .map(GroupInformation::getOneLetterCode)
                                .collect(Collectors.joining()));

                    if(!globalList.containsSequenceMotif(sequenceMotif)) {
                        globalList.addSequenceMotif(sequenceMotif);
                    }
                    logger.debug("found sequence motif: {}", sequenceMotif);
                }
            }
        });

        // set global reference to all entries
        protein.getFeatureContainer().addFeature(globalList);
    }
}