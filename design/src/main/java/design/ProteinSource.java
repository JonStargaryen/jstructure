package design;

import de.bioforscher.jstructure.feature.motif.SequenceMotif;
import de.bioforscher.jstructure.feature.motif.SequenceMotifAnnotator;
import de.bioforscher.jstructure.feature.sse.SecondaryStructureAnnotator;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import de.bioforscher.jstructure.parser.ProteinParser;
import design.parser.opm.OPMParser;
import design.parser.opm.TMHelix;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static java.util.Objects.nonNull;

/**
 * Provides access to all proteins of the data set with attached <code>*.opm</code> topology information.
 * Created by S on 31.10.2016.
 */
public class ProteinSource {
    /**
     * The universal function to load protein structures. All proteins of the data set which also provide <code>*.opm</code>
     * files will be returned with topology information attached to them.
     * @return a collection of all
     * @throws IOException when directories or files cannot be found
     */
    public static List<Protein> loadProteins() throws IOException {
        return Files.list(Paths.get(DesignConstants.OPM_RAW_DIR))
//                                      .limit(1)
                                      .map(path -> ProteinParser.parsePDBFile(DesignConstants.PDB_DIR +
                                              path.toFile().getName().split("\\.")[0] + DesignConstants.PDB_SUFFIX))
                                      .collect(Collectors.toList());
    }

    public static List<Protein> loadProteins(boolean annotateSSE, boolean annotateMotifs, boolean annotateTopology)
            throws IOException {
        List<Protein> proteins = loadProteins();

        proteins.forEach(protein -> {
            System.out.println("parsing information for " + protein.getName());
            if(annotateTopology) {
                OPMParser.parse(protein, Paths.get(DesignConstants.OPM_RAW_DIR + protein.getName() +
                        DesignConstants.OPM_SUFFIX));
            }
            if(annotateSSE) {
                addSecondaryStructureInformation(protein);
            }
            if(annotateMotifs) {
                addSequenceMotifInformation(protein);
            }
        });

        return proteins;
    }

    @SuppressWarnings("unchecked")
    public static Set<SequenceMotif> loadSequenceMotifs() throws IOException {
        Stream<SequenceMotif> motifStream = loadProteins().stream()
                .flatMap(protein -> Selection.on(protein)
                                             .aminoAcids()
                                             .asFilteredGroups())
                .filter(residue -> nonNull(residue.getFeature(List.class,
                        SequenceMotifAnnotator.FeatureNames.SEQUENCE_MOTIF)))
                .flatMap(residue -> residue.getFeature(List.class,
                        SequenceMotifAnnotator.FeatureNames.SEQUENCE_MOTIF).stream());

        return motifStream.collect(Collectors.toSet());
    }

    @SuppressWarnings("unchecked")
    public static Map<String, List<SequenceMotif>> groupSequenceMotifByTopology() throws IOException {
        // fetch all residues
        Stream<Group> residues = ProteinSource.loadProteins()
                .stream()
                .flatMap(protein -> Selection.on(protein)
                                             .aminoAcids()
                                             .asFilteredGroups());

        // create topology map
        Stream<SequenceMotif> motifs = residues
                .filter(residue -> nonNull(residue.getFeature(List.class,
                        SequenceMotifAnnotator.FeatureNames.SEQUENCE_MOTIF)))
                .flatMap(residue -> residue.getFeature(List.class,
                        SequenceMotifAnnotator.FeatureNames.SEQUENCE_MOTIF).stream());
        // for some reason mapping and instantaneously filtering will not compile
        Set<SequenceMotif> motifSet = motifs.collect(Collectors.toSet());
        // group by topology
        return motifSet.stream()
                .collect(Collectors.groupingBy(ProteinSource::determineTopologyGroup));
    }

    public static String determineTopologyGroup(SequenceMotif motif) {
        boolean startTM = nonNull(motif.getStartGroup().getFeature(TMHelix.class,
                OPMParser.FeatureNames.TM_HELIX));
        boolean endTM = nonNull(motif.getEndGroup().getFeature(TMHelix.class,
                OPMParser.FeatureNames.TM_HELIX));

        if(startTM && endTM) {
            return "I";
        }
        if(!startTM && !endTM) {
            return "o";
        }
        return "t";
    }

    private static void addSequenceMotifInformation(Protein protein) {
        SequenceMotifAnnotator sma = new SequenceMotifAnnotator();
        sma.process(protein);
    }

    private static void addSecondaryStructureInformation(Protein protein) {
        SecondaryStructureAnnotator ssa = new SecondaryStructureAnnotator();
        ssa.process(protein);
    }
}
