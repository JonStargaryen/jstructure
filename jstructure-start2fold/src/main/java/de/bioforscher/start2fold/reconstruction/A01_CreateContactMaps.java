package de.bioforscher.start2fold.reconstruction;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.feature.graphs.ProteinGraphFactory;
import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.feature.sse.SecondaryStructureType;
import de.bioforscher.jstructure.mathematics.Pair;
import de.bioforscher.jstructure.mathematics.SetOperations;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.start2fold.Start2FoldConstants;
import de.bioforscher.start2fold.model.Start2FoldResidueAnnotation;
import de.bioforscher.start2fold.parser.Start2FoldXmlParser;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

public class A01_CreateContactMaps {
    private static final Path BASE_PATH = Paths.get("/home/bittrich/git/phd_sb_repo/data/start2fold/");
    private static final Path MAP_PATH = BASE_PATH.resolve("reconstruction").resolve("maps");

    public static void main(String[] args) throws IOException {
        Files.lines(BASE_PATH.resolve("reconstruction").resolve("ids.list"))
                .forEach(A01_CreateContactMaps::handleFile);
    }

    private static void handleFile(String line) {
        String[] split = line.split(";");
        String stfId = split[0];
        String pdbId = split[1];

        List<Integer> experimentIds = Pattern.compile(",")
                .splitAsStream(split[2].replaceAll("\\[", "").replaceAll("]", ""))
                .map(Integer::valueOf)
                .collect(Collectors.toList());

        Structure structure = StructureParser.fromPath(BASE_PATH.resolve("pdb").resolve(stfId + ".pdb")).parse();
        Chain chain = structure.chainsWithAminoAcids().findFirst().get();
        String sequence = chain.getAminoAcidSequence();
        String secondaryStructureString = chain.aminoAcids()
                .map(aminoAcid -> aminoAcid.getFeature(GenericSecondaryStructure.class))
                .map(GenericSecondaryStructure::getSecondaryStructure)
                .map(SecondaryStructureType::getReducedRepresentation)
                .collect(Collectors.joining())
                .toUpperCase();

        Start2FoldConstants.write(BASE_PATH.resolve("reconstruction").resolve("fasta").resolve(stfId + ".fasta"),
                ">" + stfId + System.lineSeparator() + sequence);
        Start2FoldConstants.write(BASE_PATH.resolve("reconstruction").resolve("sse").resolve(stfId + ".sse"),
                ">" + stfId + System.lineSeparator() + secondaryStructureString);

        Start2FoldXmlParser.parseSpecificExperiment(chain,
                Start2FoldConstants.XML_DIRECTORY.resolve(stfId + ".xml"),
                experimentIds);

        List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());
        List<AminoAcid> earlyFoldingResidues = chain.aminoAcids()
                .filter(aminoAcid -> aminoAcid.getFeature(Start2FoldResidueAnnotation.class).isEarly())
                .collect(Collectors.toList());

        List<Pair<AminoAcid, AminoAcid>> contacts = SetOperations.unorderedPairsOf(aminoAcids)
                .filter(pair -> areNonCovalentGroups(pair.getLeft(), pair.getRight()))
                .filter(pair -> ProteinGraphFactory.InteractionScheme.CALPHA8.areInContact(pair.getLeft(), pair.getRight()))
                .collect(Collectors.toList());
        List<Pair<AminoAcid, AminoAcid>> earlyFoldingContacts = contacts.stream()
                .filter(pair -> earlyFoldingResidues.contains(pair.getLeft()) && earlyFoldingResidues.contains(pair.getRight()))
                .collect(Collectors.toList());

        String percentage = StandardFormat.formatToInteger(100 * earlyFoldingContacts.size() / (double) contacts.size());
        System.out.println("fraction of EFR contacts is " + percentage + "%: " + earlyFoldingContacts.size() + " " + contacts.size());

        Start2FoldConstants.write(MAP_PATH.resolve(stfId + "-sampled-100-1.rr"),
                composeRRString(contacts, sequence));
        Start2FoldConstants.write(MAP_PATH.resolve(stfId + "-efr-" + percentage + "-1.rr"),
                composeRRString(earlyFoldingContacts, sequence));

        for(int i = 5; i < 100; i = i + 5) {
            int numberOfContactsToSelect = (int) (i / (double) 100 * contacts.size());
            for(int j = 1; j < 6; j++) {
                Collections.shuffle(contacts);
                List<Pair<AminoAcid, AminoAcid>> selectedContacts = contacts.subList(0, numberOfContactsToSelect);
                Start2FoldConstants.write(MAP_PATH.resolve(stfId + "-random-" + i + "-" + j + ".rr"),
                        composeRRString(selectedContacts, sequence));
            }
        }

        // create samplings of random residues
        for(int i = 5; i < 100; i = i + 5) {
            int numberOfResiduesToSelect = (int) (i / (double) 100 * aminoAcids.size());
            for(int j = 1; j < 6; j++) {
                Collections.shuffle(aminoAcids);
                List<AminoAcid> selectedAminoAcids = aminoAcids.subList(0, numberOfResiduesToSelect);
                Start2FoldConstants.write(MAP_PATH.resolve(stfId + "-residues-" + i + "-" + j + ".rr"),
                        composeRRString(contacts.stream()
                                .filter(contact -> selectedAminoAcids.contains(contact.getLeft()) && selectedAminoAcids.contains(contact.getRight()))
                                .collect(Collectors.toList()), sequence));
            }
        }

        // create samplings of comparable nature of EFR contacts
        for(int j = 1; j < 6; j++) {
            int numberOfResiduesToSelect = earlyFoldingResidues.size();
            List<AminoAcid> interactingResidues = getInteractingResidues(aminoAcids, contacts, numberOfResiduesToSelect);
            Start2FoldConstants.write(MAP_PATH.resolve(stfId + "-interacting-" + percentage + "-" + j + ".rr"),
                    composeRRString(contacts.stream()
                            .filter(contact -> interactingResidues.contains(contact.getLeft()) && interactingResidues.contains(contact.getRight()))
                            .collect(Collectors.toList()), sequence));
        }
    }

    private static List<AminoAcid> getInteractingResidues(List<AminoAcid> aminoAcids,
                                                          List<Pair<AminoAcid, AminoAcid>> contacts,
                                                          int numberOfResiduesToSelect) {
        List<AminoAcid> selectedAminoAcids = new ArrayList<>();
        Collections.shuffle(aminoAcids);
        selectedAminoAcids.add(aminoAcids.get(0));

        while(selectedAminoAcids.size() < numberOfResiduesToSelect) {
            List<AminoAcid> adjacentAminoAcids = contacts.stream()
                    .map(contact -> mapToAdjacentAminoAcid(contact, selectedAminoAcids))
                    .filter(Optional::isPresent)
                    .map(Optional::get)
                    .collect(Collectors.toList());
            Collections.shuffle(adjacentAminoAcids);
            selectedAminoAcids.add(adjacentAminoAcids.get(0));
        }

        return selectedAminoAcids;
    }

    private static Optional<AminoAcid> mapToAdjacentAminoAcid(Pair<AminoAcid, AminoAcid> contact, List<AminoAcid> selectedAminoAcids) {
        AminoAcid aminoAcid1 = contact.getLeft();
        AminoAcid aminoAcid2 = contact.getRight();

        if(selectedAminoAcids.contains(aminoAcid1) && selectedAminoAcids.contains(aminoAcid2)) {
            return Optional.empty();
        }

        if(!selectedAminoAcids.contains(aminoAcid1) && !selectedAminoAcids.contains(aminoAcid2)) {
            return Optional.empty();
        }

        if(selectedAminoAcids.contains(aminoAcid1)) {
            return Optional.of(aminoAcid2);
        }

        return Optional.of(aminoAcid1);
    }

    private static String composeRRString(List<Pair<AminoAcid, AminoAcid>> edges, String sequence) {
        // i  j  d1  d2  p
        // 24 33 0 8 12.279515
        // i and j: residue numbers
        return edges.stream()
                .map(edge -> (edge.getLeft().getResidueIndex() + 1) + " " +
                        (edge.getRight().getResidueIndex() + 1) + " 0 8.0 1")
                .collect(Collectors.joining(System.lineSeparator(),
                        sequence + System.lineSeparator(),
                        ""));
    }

    private static boolean areNonCovalentGroups(Group group1, Group group2) {
        return Math.abs(group1.getResidueIndex() - group2.getResidueIndex()) > 1;
    }
}
