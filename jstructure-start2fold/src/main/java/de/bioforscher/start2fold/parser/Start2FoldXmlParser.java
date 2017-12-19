package de.bioforscher.start2fold.parser;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.start2fold.model.Experiment;
import de.bioforscher.start2fold.model.ProtectionLevel;
import de.bioforscher.start2fold.model.Start2FoldResidueAnnotation;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.select.Elements;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.InputStream;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class Start2FoldXmlParser {
    private static final Logger logger = LoggerFactory.getLogger(Start2FoldXmlParser.class);

    public static void parse(Chain chain, InputStream inputStream) {
        try {
            // assign baseline resp. entry container for each residue
            chain.aminoAcids().forEach(aminoAcid -> aminoAcid.getFeatureContainer().addFeature(new Start2FoldResidueAnnotation()));

            Document document = Jsoup.parse(inputStream, "UTF-8", "/");

            Elements experimentElements = document.getElementsByTag("experiment");
//            logger.info("found {} experiments",
//                    experimentElements.size());

            List<Experiment> experiments = experimentElements.stream()
                    .map(Experiment::parse)
                    .collect(Collectors.toList());

            for (Experiment experiment : experiments) {
                assignValues(experiment, chain);
            }
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    private static void assignValues(Experiment experiment, Chain chain) {
        String pdbSequence = chain.getAminoAcidSequence();
        String experimentSequence = experiment.getSequence();

        // align sequences to ensure correct mapping
        SequencePair<ProteinSequence, AminoAcidCompound> pair = null;
        try {
            pair = Alignments.getPairwiseAlignment(new ProteinSequence(experimentSequence),
                    new ProteinSequence(pdbSequence),
                    Alignments.PairwiseSequenceAlignerType.GLOBAL,
                    new SimpleGapPenalty(),
                    SubstitutionMatrixHelper.getBlosum62());

//            logger.info("alignment:{}{}",
//                    System.lineSeparator(),
//                    pair.toString());

            List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());

            for(Experiment.Residue residue : experiment.getResidues()) {
                int experimentIndex = residue.getIndex() - 1;
                try {
                    int pdbIndex;
                    if(residue.getCode().equals("P") && residue.getIndex() == 1) {
                        // super-russian fix for STF0017 where the alignment should match theoretically
                        pdbIndex = 0;
                    } else {
                        pdbIndex = pair.getIndexInTargetForQueryAt(experimentIndex);
                    }
                    AminoAcid aminoAcid = aminoAcids.get(pdbIndex);

//                logger.debug("mapped experiment {}-{} onto PDB {}-{}",
//                        residue.getCode(),
//                        residue.getIndex(),
//                        aminoAcid.getOneLetterCode(),
//                        aminoAcid.getResidueIdentifier());

//                if(!residue.getCode().equals(aminoAcid.getOneLetterCode())) {
//                    logger.warn("alignment:{}{}",
//                            System.lineSeparator(),
//                            pair.toString());
//                    logger.warn("could not map data correctly for {}-{} and {}-{}",
//                            residue.getCode(),
//                            residue.getIndex(),
//                            aminoAcid.getOneLetterCode(),
//                            aminoAcid.getResidueIdentifier());
//                }

                    // ignore: STF0034 (cannot align)

                    // assign experiment-specific protection level to residue
                    aminoAcid.getFeature(Start2FoldResidueAnnotation.class).addProtectionLevelEntry(Stream.of(ProtectionLevel.values())
                            .filter(pl -> pl == experiment.getProtectionLevel())
                            .findFirst()
                            .get());
                } catch (Exception e) {
                    // residue not present in structure - e.g. for STF0031 and STF0032
                    logger.warn("alignment:{}{}",
                            System.lineSeparator(),
                            pair.toString());
                    logger.warn("failed to map residue {}-{}",
                            residue.getCode(),
                            residue.getIndex(),
                            e);
                }
            }
        } catch (CompoundNotFoundException e) {
            throw new IllegalArgumentException(e);
        }
    }

    private static void assignValuesForEarly(Experiment experiment, Chain chain) {
        String pdbSequence = chain.getAminoAcidSequence();
        String experimentSequence = experiment.getSequence();

        // align sequences to ensure correct mapping
        SequencePair<ProteinSequence, AminoAcidCompound> pair = null;
        try {
            pair = Alignments.getPairwiseAlignment(new ProteinSequence(experimentSequence),
                    new ProteinSequence(pdbSequence),
                    Alignments.PairwiseSequenceAlignerType.GLOBAL,
                    new SimpleGapPenalty(),
                    SubstitutionMatrixHelper.getBlosum62());

            List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());

            for(Experiment.Residue residue : experiment.getResidues()) {
                int experimentIndex = residue.getIndex() - 1;
                try {
                    int pdbIndex;
                    if(residue.getCode().equals("P") && residue.getIndex() == 1) {
                        // super-russian fix for STF0017 where the alignment should match theoretically
                        pdbIndex = 0;
                    } else {
                        pdbIndex = pair.getIndexInTargetForQueryAt(experimentIndex);
                    }
                    AminoAcid aminoAcid = aminoAcids.get(pdbIndex);

                    // assign experiment-specific protection level to residue
                    aminoAcid.getFeature(Start2FoldResidueAnnotation.class).addProtectionLevelEntry(ProtectionLevel.EARLY);
                } catch (Exception e) {
                    // residue not present in structure - e.g. for STF0031 and STF0032
                    logger.warn("alignment:{}{}",
                            System.lineSeparator(),
                            pair.toString());
                    logger.warn("failed to map residue {}-{}",
                            residue.getCode(),
                            residue.getIndex(),
                            e);
                }
            }
        } catch (CompoundNotFoundException e) {
            throw new IllegalArgumentException(e);
        }
    }

    public static void parse(Chain chain, Path path) {
        try {
            parse(chain, Files.newInputStream(path));
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    public static void parseSpecificExperiment(Chain chain, InputStream inputStream, List<Integer> experimentIds) {
//        logger.info("ids: {}",
//                experimentIds);
        try {
            // assign baseline resp. entry container for each residue
            chain.aminoAcids().forEach(aminoAcid -> aminoAcid.getFeatureContainer().addFeature(new Start2FoldResidueAnnotation()));

            Document document = Jsoup.parse(inputStream, "UTF-8", "/");

            Elements experimentElements = document.getElementsByTag("experiment");

            List<Experiment> experiments = experimentElements.stream()
                    .map(Experiment::parse)
                    .filter(experiment -> experimentIds.contains(experiment.getExperimentId()))
                    .collect(Collectors.toList());

//            logger.info("processing {} experiments", experiments.size());

            if(experimentIds.size() != experiments.size()) {
                throw new IllegalArgumentException("did not find all referenced experiments");
            }

            for (Experiment experiment : experiments) {
                assignValuesForEarly(experiment, chain);
            }
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    public static void parseSpecificExperiment(Chain chain, Path path, List<Integer> experimentIds) {
        try {
            parseSpecificExperiment(chain, Files.newInputStream(path), experimentIds);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }
}
