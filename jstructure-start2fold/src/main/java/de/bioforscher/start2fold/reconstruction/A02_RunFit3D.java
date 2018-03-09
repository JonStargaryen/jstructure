package de.bioforscher.start2fold.reconstruction;

import de.bioforscher.singa.structure.algorithms.superimposition.fit3d.Fit3D;
import de.bioforscher.singa.structure.algorithms.superimposition.fit3d.Fit3DBuilder;
import de.bioforscher.singa.structure.model.interfaces.Structure;
import de.bioforscher.singa.structure.model.oak.StructuralEntityFilter;
import de.bioforscher.singa.structure.model.oak.StructuralMotif;
import de.bioforscher.singa.structure.parser.pdb.structures.StructureParser;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

public class A02_RunFit3D {
    private static final Path BASE_PATH = Paths.get("/home/bittrich/git/phd_sb_repo/data/start2fold/reconstruction/fit3d/");

    public static void main(String[] args) throws IOException {
        Files.list(BASE_PATH)
                .filter(path -> path.toFile().getName().contains("target"))
                .filter(path -> path.toFile().getName().contains("2lzm"))
                .forEach(A02_RunFit3D::handleTarget);
    }

    private static void handleTarget(Path target) {
        try {
            String pdbId = target.toFile().getName().split("-")[0];
            Path motif = BASE_PATH.resolve(pdbId + "-motif.pdb");
            System.out.println("handling " + target + " vs " + motif);

            Structure motifStructure = StructureParser.local().path(motif).parse();
            Structure targetStructure = StructureParser.local().path(target).parse();

            StructuralMotif queryMotif = StructuralMotif.fromLeafSubstructures(motifStructure.getAllLeafSubstructures());

            Fit3D run = Fit3DBuilder.create()
                    .query(queryMotif)
                    .target(targetStructure)
                    .atomFilter(StructuralEntityFilter.AtomFilter.isArbitrary())
                    .rmsdCutoff(2.0)
                    .run();

            run.writeSummaryFile(Paths.get("/home/bittrich/fit3d/output.fit3d"));
            run.writeMatches(Paths.get("/home/bittrich/fit3d/matches"), .75);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }
}
