package de.bioforscher.start2fold.reconstruction;

import de.bioforscher.jstructure.feature.graphs.ProteinGraphFactory;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.mathematics.Pair;
import de.bioforscher.jstructure.mathematics.SetOperations;
import de.bioforscher.jstructure.model.identifier.IdentifierFactory;
import de.bioforscher.jstructure.model.structure.*;
import de.bioforscher.jstructure.model.structure.aminoacid.Alanine;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.start2fold.Start2FoldConstants;
import de.bioforscher.start2fold.model.Start2FoldResidueAnnotation;
import de.bioforscher.start2fold.parser.Start2FoldXmlParser;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

public class A01_WriteContactRepresentationFiles {
    private static final Logger logger = LoggerFactory.getLogger(A01_WriteContactRepresentationFiles.class);
    private static int counter;

    public static void main(String[] args) throws IOException {
        Files.lines(Start2FoldConstants.PANCSA_LIST)
                .forEach(A01_WriteContactRepresentationFiles::handleLine);
    }

    private static void handleLine(String line) {
        try {
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

            List<AminoAcid> earlyFoldingResidues = chain.aminoAcids()
                    .filter(aminoAcid -> aminoAcid.getFeature(Start2FoldResidueAnnotation.class).isEarly())
                    .collect(Collectors.toList());

            // filter for contacts of EFR
            PLIPInteractionContainer plipInteractionContainer = chain.getFeature(PLIPInteractionContainer.class);

            // ignore really small or non-interacting groups of EFR
            if(plipInteractionContainer.getInteractions()
                    .stream()
                    .filter(interaction -> earlyFoldingResidues.contains(interaction.getPartner1()) && earlyFoldingResidues.contains(interaction.getPartner2()))
                    .count() < 2) {
                return;
            }

            System.out.println(line);

            List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());

            counter = 1;
            String efrOutput = structure.getHeader() + createPdbRepresentation(earlyFoldingResidues);
            counter = 1;
            String fullOutput = structure.getHeader() + createPdbRepresentation(aminoAcids);

            Files.write(Paths.get("/home/bittrich/git/phd_sb_repo/data/start2fold/reconstruction/fit3d/" +
                    pdbId + "-motif.pdb"),
                    efrOutput.getBytes());
            Files.write(Paths.get("/home/bittrich/git/phd_sb_repo/data/start2fold/reconstruction/fit3d/" +
                            pdbId + "-target.pdb"),
                    fullOutput.getBytes());
        } catch (Exception e) {
            e.printStackTrace();
            logger.info("calculation failed for {}",
                    line,
                    e);
        }
    }

    private static String createPdbRepresentation(List<AminoAcid> aminoAcids) {
        return SetOperations.unorderedPairsOf(aminoAcids)
                .filter(pair -> ProteinGraphFactory.InteractionScheme.CALPHA8.areInContact(pair.getLeft(), pair.getRight()))
                .map(A01_WriteContactRepresentationFiles::createAtom)
                .map(Atom::getPdbRepresentation)
                .collect(Collectors.joining(System.lineSeparator()));
    }

    private static Atom createAtom(Pair<AminoAcid, AminoAcid> pair) {
        double[] coordinates1 = pair.getLeft().getCa().getCoordinates();
        double[] coordinates2 = pair.getRight().getCa().getCoordinates();

        Atom atom = Atom.builder(Element.C, new double[] {
                0.5 * (coordinates1[0] + coordinates2[0]),
                0.5 * (coordinates1[1] + coordinates2[1]),
                0.5 * (coordinates1[2] + coordinates2[2])
        }).pdbSerial(counter).name("CA").build();

        Group group = new Alanine(IdentifierFactory.createResidueIdentifier(counter));
        group.addAtom(atom);

        counter++;
        return atom;
    }
}
