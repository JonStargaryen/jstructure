package de.bioforscher.jstructure.contacts;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.testutil.FileUtils;

import java.nio.file.Path;
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
