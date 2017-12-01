package de.bioforscher.jstructure.contacts;

import de.bioforscher.jstructure.mathematics.Pair;
import de.bioforscher.jstructure.mathematics.SetOperations;
import de.bioforscher.jstructure.mathematics.graph.Edge;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.testutil.FileUtils;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class ContactsConstants extends FileUtils {
    public static final Path RECONSTRUCTION_DIRECTORY = DATA_DIRECTORY.resolve("reconstruction");
    public static final Path START2FOLD_DIRECTORY = DATA_DIRECTORY.resolve("reconstruction-start2fold");

    public static Atom getBetaCarbon(AminoAcid aminoAcid) {
        return aminoAcid.atoms()
                .filter(atom -> atom.getName().equals("CB"))
                .findFirst()
                .orElse(aminoAcid.getCa());
    }

    public static TMAlignResult parseTMAlignResultFile(Path tmAlignPath) {
        return new TMAlignResult(tmAlignPath);
    }

    private static final double BETA_THRESHOLD = 8.0;
    public static List<Edge<AminoAcid>> determineNaiveInteractions(List<AminoAcid> aminoAcids) {
        List<Edge<AminoAcid>> naiveEdges = new ArrayList<>();
        List<Pair<AminoAcid, AminoAcid>> pairs = SetOperations.uniquePairsOf(aminoAcids).collect(Collectors.toList());

        for(Pair<AminoAcid, AminoAcid> pair : pairs) {
            AminoAcid aminoAcid1 = pair.getLeft();
            AminoAcid aminoAcid2 = pair.getRight();

            if(!areNonCovalentGroups(aminoAcid1, aminoAcid2)) {
                continue;
            }

            // conventional interaction criterion
            if(ContactsConstants.getBetaCarbon(aminoAcid1).calculate().distance(ContactsConstants.getBetaCarbon(aminoAcid2)) < BETA_THRESHOLD) {
                naiveEdges.add(new Edge<>(aminoAcid1, aminoAcid2, BETA_THRESHOLD));
            }
        }

        return naiveEdges;
    }

    public static boolean areNonCovalentGroups(Group group1, Group group2) {
        return Math.abs(group1.getResidueIndex() - group2.getResidueIndex()) > 1;
    }

    public static class TMAlignResult {
        private final double rmsd;
        private final double tmscore;

        TMAlignResult(Path tmAlignPath) {
            List<String> lines;
            try(Stream<String> stream = ContactsConstants.lines(tmAlignPath)) {
                lines = stream.collect(Collectors.toList());
            }
            double rmsd = -1;
            double tmscore = -1;
            for(String line : lines) {
                if(line.startsWith("Aligned length=")) {
                    rmsd = Double.valueOf(line.split("RMSD=")[1].trim().split(",")[0].trim());
                }
                if(line.startsWith("TM-score=")) {
                    tmscore = Double.valueOf(line.split("TM-score=")[1].trim().split("\\(")[0].trim());
                }
            }

            if(rmsd == -1 || tmscore == -1) {
                throw new IllegalArgumentException("tmalign result file malformed - " + tmAlignPath.toFile().getAbsolutePath());
            }

            this.rmsd = rmsd;
            this.tmscore = tmscore;
        }

        public double getRmsd() {
            return rmsd;
        }

        public double getTmscore() {
            return tmscore;
        }
    }
}
