package de.bioforscher.start2fold.parser;

import de.bioforscher.jstructure.feature.mapping.ResidueMapping;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.start2fold.model.FunctionalResidueAnnotation;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.List;

public class FunctionalResidueParser {
    private static final Logger logger = LoggerFactory.getLogger(FunctionalResidueParser.class);

//    public static void parse(Chain chain, String uniProtId, List<Integer> functionalResidues) {
//        try {
//            System.out.println(uniProtId);
//            // assign baseline resp. entry container for each residue
//            chain.aminoAcids().forEach(aminoAcid -> aminoAcid.getFeatureContainer().addFeature(new FunctionalResidueAnnotation()));
//
//            Document document = Jsoup.connect("https://www.uniprot.org/uniprot/" + uniProtId + ".xml").get();
//
//            String uniProtSequence = document.getElementsByTag("sequence")
//                    .first()
//                    .text()
//                    .replaceAll("\\s+", "");
//            String pdbSequence = chain.getAminoAcidSequence();
//
//            System.out.println(uniProtSequence);
//            System.out.println(pdbSequence);
//
//            // align sequences to ensure correct mapping
//            SequencePair<ProteinSequence, AminoAcidCompound> pair;
//            try {
//                pair = Alignments.getPairwiseAlignment(new ProteinSequence(uniProtSequence),
//                        new ProteinSequence(pdbSequence),
//                        Alignments.PairwiseSequenceAlignerType.GLOBAL,
//                        new SimpleGapPenalty(),
//                        SubstitutionMatrixHelper.getBlosum62());
//
//                List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());
//
//                for(int functionalResidue : functionalResidues) {
//                    try {
//                        int pdbIndex = pair.getIndexInTargetForQueryAt(functionalResidue);
//                        AminoAcid aminoAcid = aminoAcids.get(pdbIndex);
//
//                        aminoAcid.getFeature(FunctionalResidueAnnotation.class).addFunctionalAnnotation("functional");
//                    } catch (Exception e) {
//                        logger.warn("alignment:{}{}",
//                                System.lineSeparator(),
//                                pair.toString());
//                    }
//                }
//            } catch (CompoundNotFoundException e) {
//                throw new IllegalArgumentException(e);
//            }
//        } catch (IOException e) {
//            throw new UncheckedIOException(e);
//        }
//    }

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
