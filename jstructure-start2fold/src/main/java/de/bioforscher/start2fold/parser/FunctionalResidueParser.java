package de.bioforscher.start2fold.parser;

import de.bioforscher.jstructure.feature.mapping.ResidueMapping;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.start2fold.model.FunctionalResidueAnnotation;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.List;

public class FunctionalResidueParser {
    private static final Logger logger = LoggerFactory.getLogger(FunctionalResidueParser.class);

    public static void parse(Chain chain, List<Integer> functionalResidues) {
        // assign baseline resp. entry container for each residue
        chain.aminoAcids().forEach(aminoAcid -> aminoAcid.getFeatureContainer().addFeature(new FunctionalResidueAnnotation()));

        for(int functionalResidue : functionalResidues) {
            chain.aminoAcids()
                    .filter(aminoAcid -> aminoAcid.getFeature(ResidueMapping.class)
                            .getUniProtResidueNumber()
                            .equals(String.valueOf(functionalResidue)))
                    .forEach(aminoAcid -> aminoAcid.getFeature(FunctionalResidueAnnotation.class)
                            .addFunctionalAnnotation("functional"));
        }
    }
}
