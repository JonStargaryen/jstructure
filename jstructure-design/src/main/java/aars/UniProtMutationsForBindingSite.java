package aars;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.ResidueNumber;
import de.bioforscher.jstructure.model.structure.identifier.ProteinIdentifier;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * Search for mutations for all residues in the extracted binding site.
 * Created by bittrich on 3/28/17.
 */
class UniProtMutationsForBindingSite {
    private static final Path EFFECTS_TSV = Paths.get("/home/bittrich/git/aars_analysis/data/mutations/mutation-effects.tsv");
    private static final StringBuilder output = new StringBuilder();
    private static final StringBuilder warnings = new StringBuilder();
    private static final String DEL = "\t";
    private static final String classToProcess = "C2";

    public static void main(String[] args) throws IOException {
        output.append("pdbId" + DEL + "chainId" + DEL + "uniProtId" + DEL + "class" + DEL + "aa" + DEL + "mode" + DEL + "renumPos" + DEL + "origPos" + DEL + "type" + DEL + "uniProtPos" + DEL + "original" + DEL + "variant" + DEL + "description").append(System.lineSeparator());
        Files.walk(Paths.get("/home/bittrich/git/aars_analysis/data/binding_sites/" + classToProcess + "/ligand_based/per_type/"))
                // ignore additional files
                .filter(path -> path.toFile().isDirectory())
                // only consider truncated/curated files
                .filter(path -> path.toFile().getName().contains("truncated"))
                .flatMap(AARSConstants::list)
                .map(path -> ProteinParser.source(path).forceProteinName(ProteinIdentifier.createFromAdditionalName(path.toFile().getAbsolutePath())).parse())
                .map(UniProtMutationsForBindingSite::extractBindingSiteResidueNumbers)
                .filter(UniProtMutationsForBindingSite::hasMutations)
                .forEach(UniProtMutationsForBindingSite::handleBindingSite);

        Files.write(Paths.get("/home/bittrich/git/aars_analysis/data/mutations/binding-site-effects-" + classToProcess + ".tsv"), output.toString().getBytes());
        System.err.println(warnings.toString());
    }

    private static boolean hasMutations(BindingSite bindingSite) {
        return AARSConstants.lines(EFFECTS_TSV)
                .map(line -> line.split("\t"))
                .filter(split -> split[0].equals(bindingSite.pdbId))
                .anyMatch(split -> split[1].equals(bindingSite.chainId));
    }

    private static void handleBindingSite(BindingSite bindingSite) {
        String uniProtId = AARSConstants.lines(EFFECTS_TSV)
                .map(line -> line.split(DEL))
                .filter(split -> split[0].equals(bindingSite.pdbId))
                .filter(split -> split[1].equals(bindingSite.chainId))
                .findAny()
                .get()[2];

        // load original, full structure
        Chain originalChain = ProteinParser.source(bindingSite.pdbId).parse()
                .select()
                .chainName(bindingSite.chainId)
                .asChain();
        String pdbSequence = originalChain.getAminoAcidSequence();
        String uniProtSequence = loadUniProtSequence(uniProtId);

        // align sequences
        SequencePair<ProteinSequence, AminoAcidCompound> alignment = needle(uniProtSequence, pdbSequence);

        System.out.println(bindingSite);
        System.out.println(alignment);

        // load renumbered, but not transformed chain
        Chain renumberedChain = ProteinParser.source(Paths.get("/home/bittrich/git/aars_analysis/data/msa/" + classToProcess + "/renumbered_structures/" + bindingSite.pdbId + "_renum.pdb")).parse()
                .select()
                .chainName(bindingSite.chainId)
                .asChain();

        // key: renumbered, transformed binding site group - value: original group in PDB chain
        List<Integer> residueIndices = bindingSite.residues.stream()
                .map(Group::getResidueNumber)
                .map(ResidueNumber::getResidueNumber)
                .collect(Collectors.toList());
        Map<Group, Group> groupMapping = renumberedChain.aminoAcids()
                .filter(aminoAcid -> residueIndices.contains(aminoAcid.getResidueNumber()))
                .collect(Collectors.toMap(Function.identity(),
                // map each group to the entity in the not renumbered structure
                renumberedGroup -> originalChain.select()
                        .groupName(renumberedGroup.getThreeLetterCode())
                        .asFilteredGroups()
                        .min(Comparator.comparingDouble(originalGroup -> originalGroup.calculate().centroid().distanceFast(renumberedGroup.calculate().centroid())))
                        .get()));

        // determine sequence position in sequence alignment - rarely these indices do not match
        groupMapping.entrySet().forEach(entry -> {
            // map to index in pdb entry
            int residueIndex = originalChain.getGroups().indexOf(entry.getValue()) + 1;

            // map to index in alignment - query is the uniProtSequence - however, getIndexInTarget seems to return the correct position
            System.out.print("mapped: " + entry.getKey().getIdentifier() + " -> " +
                    entry.getValue().getIdentifier() + " -> uniprot ");
            try {
                int indexInUniProt = alignment.getIndexInQueryForTargetAt(residueIndex);
                System.out.println(alignment.getCompoundInQueryAt(indexInUniProt).getLongName().toUpperCase() + "-" + indexInUniProt);

                String indexToFind = String.valueOf(indexInUniProt);
                AARSConstants.lines(EFFECTS_TSV)
                        .map(line -> line.split("\t"))
                        .filter(split -> split[0].equals(bindingSite.pdbId))
                        .filter(split -> split[1].equals(bindingSite.chainId))
                        .filter(split -> refersToPosition(split, indexToFind))
                        .forEach(split -> {
                            // compose output line
                            String outputLine = bindingSite.pdbId + DEL +
                                    bindingSite.chainId + DEL +
                                    split[2] + DEL +
                                    bindingSite.clazz + DEL +
                                    bindingSite.aa + DEL +
                                    bindingSite.mode + DEL +
                                    entry.getKey().getResidueNumber() + DEL +
                                    entry.getValue().getResidueNumber() + DEL +
                                    split[3] + DEL +
                                    split[4] + DEL +
                                    split[5] + DEL +
                                    split[6] + DEL +
                                    split[7] + System.lineSeparator();
                            System.out.println(outputLine);
                            output.append(outputLine);
                        });
            } catch (ArrayIndexOutOfBoundsException e) {
                System.out.println("failed!");
                warnings.append("#could not map ").append(entry.getValue().getIdentifier()).append(" in ").append(bindingSite.pdbId).append("_").append(bindingSite.chainId).append(" to UniProt sequence").append(System.lineSeparator());
            }
        });
    }

    private static final Pattern positionPattern = Pattern.compile(",");

    private static boolean refersToPosition(String[] split, String indexToFind) {
        return positionPattern.splitAsStream(split[4].replace("[", "").replace("]", ""))
                .map(String::trim)
                .anyMatch(position -> position.equals(indexToFind));
    }

    private static String loadUniProtSequence(String uniProtId) {
        try {
            Document document = Jsoup.connect("http://www.uniprot.org/uniprot/" + uniProtId + ".xml").get();
            return document.getElementsByTag("sequence").text().replaceAll("\\s+", "");
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    private static BindingSite extractBindingSiteResidueNumbers(Protein protein) {
        return new BindingSite(protein.getPdbId().getPdbId(),
                protein.chains().findFirst().get().getChainId().getChainId(),
                protein.getGroups());
    }

    private static SequencePair<ProteinSequence, AminoAcidCompound> needle(String sequence1, String sequence2) {
        try {
            return Alignments.getPairwiseAlignment(new ProteinSequence(sequence1),
                    new ProteinSequence(sequence2),
                    Alignments.PairwiseSequenceAlignerType.GLOBAL,
                    new SimpleGapPenalty(),
                    SubstitutionMatrixHelper.getBlosum62());
        } catch (CompoundNotFoundException e) {
            throw new IllegalArgumentException(e);
        }
    }

    static class BindingSite {
        String pdbId, chainId, source;
        String clazz, aa, mode;
        List<Group> residues;

        BindingSite(String source, String chainId, List<Group> residues) {
            File file = new File(source);
            // name is something cryptic like: truncated_3mf1_A_renum
            this.pdbId = file.getName().split("_")[1];
            this.chainId = chainId;
            this.source = source;
            this.residues = residues;
            String[] split = file.getAbsolutePath().split("/");
            this.clazz = split[7];
            this.aa = split[11];
            this.mode = split[10];
        }

        @Override
        public String toString() {
            return "pdbId=" + pdbId + " chainId=" + chainId + " source=" + source + " residues=" + residues.stream()
                    .map(Group::getResidueNumber)
                    .collect(Collectors.toList());
        }
    }
}
