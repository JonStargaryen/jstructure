package de.bioforscher.start2fold.collect;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.feature.interactions.PLIPInteraction;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.start2fold.Start2FoldConstants;
import de.bioforscher.start2fold.model.Start2FoldResidueAnnotation;
import de.bioforscher.start2fold.parser.Start2FoldXmlParser;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Optional;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * Characterize and count contacts between EFR with other EFR. Can those be the basis of reconstructions?
 */
public class A05_WriteEarlyFoldingContactCsv {
    private static final Logger logger = LoggerFactory.getLogger(A05_WriteEarlyFoldingContactCsv.class);
    private static List<AminoAcid> earlyFoldingResidues;

    public static void main(String[] args) throws IOException {
        String localOutput = Files.lines(Start2FoldConstants.PANCSA_LIST)
                .map(A05_WriteEarlyFoldingContactCsv::handleLineLocally)
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.joining(System.lineSeparator(),
                        "pdb,chain,aa1,aa2,sse1,sse2," +
                                "type,distance,efr" + System.lineSeparator(),
                        ""));

        Start2FoldConstants.write(Start2FoldConstants.STATISTICS_DIRECTORY.resolve("efr-contacts.csv"),
                localOutput);
    }

    private static Optional<String> handleLineLocally(String line) {
        try {
            System.out.println(line);
            String[] split = line.split(";");
            String entryId = split[0];
            String pdbId = split[1];
            List<Integer> experimentIds = Pattern.compile(",")
                    .splitAsStream(split[2].replaceAll("\\[", "").replaceAll("]", ""))
                    .map(Integer::valueOf)
                    .collect(Collectors.toList());

            Structure structure = StructureParser.fromPdbId(pdbId).parse();
            Chain chain = structure.chains().findFirst().get();

            Start2FoldXmlParser.parseSpecificExperiment(chain,
                    Start2FoldConstants.XML_DIRECTORY.resolve(entryId + ".xml"),
                    experimentIds);

            earlyFoldingResidues = chain.aminoAcids()
                    .filter(aminoAcid -> aminoAcid.getFeature(Start2FoldResidueAnnotation.class).isEarly())
                    .collect(Collectors.toList());

            // filter for contacts of EFR
            PLIPInteractionContainer plipInteractionContainer = chain.getFeature(PLIPInteractionContainer.class);

            System.out.println(pdbId + "\t" +
                    earlyFoldingResidues.size() + "\t" +
                    StandardFormat.format(plipInteractionContainer.getInteractions().size() / (double) earlyFoldingResidues.size()) + "\t" +
                    plipInteractionContainer.getInteractions().size());

            // ignore really small or non-interacting groups of EFR
            if(plipInteractionContainer.getInteractions().size() < 2) {
                return Optional.empty();
            }

            return Optional.of(plipInteractionContainer.getInteractions()
                    .stream()
                    .map(A05_WriteEarlyFoldingContactCsv::handleInteraction)
                    .collect(Collectors.joining(System.lineSeparator())));
        } catch (Exception e) {
            e.printStackTrace();
            logger.info("calculation failed for {}",
                    line,
                    e);
            return Optional.empty();
        }
    }

    private static String handleInteraction(PLIPInteraction plipInteraction) {
        AminoAcid aa1 = (AminoAcid) plipInteraction.getPartner1();
        AminoAcid aa2 = (AminoAcid) plipInteraction.getPartner2();
        String pdbId = aa1.getParentStructure().getProteinIdentifier().getPdbId();
        return pdbId + "," +
                "A" + "," +
                aa1.getOneLetterCode() + "," +
                aa1.getFeature(GenericSecondaryStructure.class).getSecondaryStructure().getOneLetterRepresentation() + "," +
                aa2.getOneLetterCode() + "," +
                aa2.getFeature(GenericSecondaryStructure.class).getSecondaryStructure().getOneLetterRepresentation() + "," +
                plipInteraction.getClass().getSimpleName() + "," +
                (Math.abs(aa1.getResidueIndex() - aa2.getResidueIndex()) > 5 ? "long-range" : "local") + "," +
                (earlyFoldingResidues.contains(aa1) && earlyFoldingResidues.contains(aa2) ? "efr" : "normal");
    }
}
