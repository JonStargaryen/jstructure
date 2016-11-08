package design.aggregation;

import de.bioforscher.jstructure.alignment.AlignmentAlgorithm;
import de.bioforscher.jstructure.alignment.AlignmentResult;
import de.bioforscher.jstructure.alignment.SVDSuperimposer;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.filter.AtomNameFilter;
import de.bioforscher.jstructure.parser.ProteinParser;
import design.DesignConstants;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;

/**
 * Aligns fragments of one sequence motif observation.
 * Created by S on 07.11.2016.
 */
public class S07_AlignFragments {
    public static void main(String[] args) {
        S06_ExtractSequences.TOPOLOGIES.stream()
                .limit(1)
                .forEach(S07_AlignFragments::handleTopology);
    }

    static class ProteinSuperimposerByReference implements Consumer<Protein> {
        /**
         * the protein all entries are aligned against
         */
        Protein reference;
        /**
         * the container of aligned proteins
         */
        List<Protein> alignedProteins;
        /**
         * the approach to align fragments
         */
        AlignmentAlgorithm alignmentStrategy;

        ProteinSuperimposerByReference() {
            alignedProteins = new ArrayList<>();
            alignmentStrategy = new SVDSuperimposer(AtomNameFilter.BACKBONE_ATOM_FILTER);
        }

        @Override
        public void accept(Protein protein) {
            if(reference == null) {
                reference = protein;
            }

            alignedProteins.add(protein);
            AlignmentResult alignmentResult = alignmentStrategy.align(reference, protein);
            System.out.printf("%2f A is the rmsd between %s and %s%n", alignmentResult.getRmsd(),
                    reference.getName(), protein.getName());
            alignmentResult.transform(protein);
        }

        void combine(ProteinSuperimposerByReference other) {
            other.alignedProteins.forEach(this);
        }

        List<Protein> getAlignedProteins() {
            return alignedProteins;
        }
    }

    private static void handleTopology(final String topology) {
        final String motif = "AL6";
        try {
            List<Protein> alignedProteins = Files.list(Paths.get(DesignConstants.MOTIF_FRAGMENT_BY_TOPOLOGY_DIR + topology + "/"))
                .filter(path -> path.getFileName().toString().startsWith(motif))
                .map(ProteinParser::parsePDBFile)
//                .limit(5)
                .collect(ProteinSuperimposerByReference::new,
                      ProteinSuperimposerByReference::accept,
                      ProteinSuperimposerByReference::combine)
                .getAlignedProteins();
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }
}
