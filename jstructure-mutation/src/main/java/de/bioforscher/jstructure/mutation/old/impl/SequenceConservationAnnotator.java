package de.bioforscher.jstructure.mutation.old.impl;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.util.Map;
import java.util.Optional;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * A custom feature provider which annotates the sequence conservation profile for each amino acid.
 * Created by bittrich on 7/12/17.
 */
@Deprecated
class SequenceConservationAnnotator {
    static void process(Chain chain, Map<String, String> alignmentMap) {
        int alignmentLength = alignmentMap.values().stream().findFirst().get().length();
        double frequencyIncrease = 1 / (double) alignmentMap.size();

        for(int position = 0; position < alignmentLength; position++) {
            // first: check whether group is present in chain
            Optional<AminoAcid> aminoAcidOptional = chain.select()
                    .residueNumber(position + 1)
                    .asOptionalAminoAcid();
            if(!aminoAcidOptional.isPresent()) {
                continue;
            }

            double deletionFrequency = 0.0;
            Map<AminoAcid.Family, Double> frequencyTable = Stream.of(AminoAcid.Family.values())
                    .collect(Collectors.toMap(Function.identity(), aminoAcid -> 0.0));
            for(String alignment : alignmentMap.values()) {
                char characterInAlignment = alignment.charAt(position);

                // count deletions
                if(characterInAlignment == '-') {
                    deletionFrequency += frequencyIncrease;
                    continue;
                }

                // increase frequency of normal amino acids
                AminoAcid.Family aminoAcid = AminoAcid.Family.resolveOneLetterCode(characterInAlignment);
                //TODO improve by natural distribution of amino acids
                frequencyTable.put(aminoAcid, frequencyTable.get(aminoAcid) + frequencyIncrease);
            }

            SequenceConservationProfile sequenceConservationProfile = new SequenceConservationProfile(null, frequencyTable, deletionFrequency);
            aminoAcidOptional.get().getFeatureContainer().addFeature(sequenceConservationProfile);
        }
    }
}
