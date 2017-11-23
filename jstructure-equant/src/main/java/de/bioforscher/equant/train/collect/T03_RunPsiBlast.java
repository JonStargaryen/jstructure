package de.bioforscher.equant.train.collect;

import de.bioforscher.equant.EquantConstants;
import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.align.impl.LocalBlastWrapper;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.nio.file.Path;
import java.util.stream.Collectors;

/**
 * For all targets run PSI-BLAST and generate evolutionary information profiles.
 */
public class T03_RunPsiBlast {
    private static final Logger logger = LoggerFactory.getLogger(T03_RunPsiBlast.class);
    private static final Path DIRECTORY = EquantConstants.CASP9_DIRECTORY;
    private static final Path OUTPUT_DIRECTORY = DIRECTORY.resolve("psiblast");
    private static final LocalBlastWrapper LOCAL_BLAST_WRAPPER = new LocalBlastWrapper();

    public static void main(String[] args) {
        EquantConstants.list(DIRECTORY.resolve("fasta"))
                .forEach(T03_RunPsiBlast::createPsiBlastProfile);
    }

    private static void createPsiBlastProfile(Path path) {
        String name = path.toFile().getName().split("\\.")[0];
        String sequence = EquantConstants.lines(path)
                .filter(line -> !line.startsWith(">"))
                .collect(Collectors.joining());
        logger.info("handling {} with sequence:{}{}",
                name,
                System.lineSeparator(),
                sequence);
        LocalBlastWrapper.PsiBlastResult psiBlastResult = LOCAL_BLAST_WRAPPER.executePsiBlastUniref50(sequence);
        String psiBlastProfile = psiBlastResult.getInformation().stream()
                .map(StandardFormat::format)
                .collect(Collectors.joining(System.lineSeparator()));
        EquantConstants.write(OUTPUT_DIRECTORY.resolve(name + ".psiblast"), psiBlastProfile);
    }
}
