package de.bioforscher.explorer.membrane;

import de.bioforscher.explorer.membrane.model.representative.ExplorerProtein;
import de.bioforscher.jstructure.alignment.SVDSuperimposer;
import de.bioforscher.jstructure.alignment.StructureAlignmentResult;
import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceAreaCalculator;
import de.bioforscher.jstructure.feature.cerosene.SequenceCerosene;
import de.bioforscher.jstructure.feature.motif.SequenceMotifAnnotator;
import de.bioforscher.jstructure.feature.sse.SecondaryStructureAnnotator;
import de.bioforscher.jstructure.feature.topology.ANVIL;
import de.bioforscher.jstructure.feature.topology.Membrane;
import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProviderRegistry;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.parser.ProteinParser;
import de.bioforscher.jstructure.parser.opm.OPMDatabaseQuery;
import de.bioforscher.jstructure.parser.plip.PLIPAnnotator;
import de.bioforscher.jstructure.parser.sifts.SiftsParser;
import de.bioforscher.jstructure.parser.uniprot.UniProtAnnotator;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import javax.annotation.PostConstruct;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Objects;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Access to the database service.
 * Created by bittrich on 2/20/17.
 */
@Service
public class ProteinService {
    private static final Logger logger = LoggerFactory.getLogger(ProteinService.class);
    private ProteinRepository repository;
    private List<String> allProteins;
    private List<AbstractFeatureProvider> featureProviders = Stream.of(AccessibleSurfaceAreaCalculator.RELATIVE_ACCESSIBLE_SURFACE_AREA,
            SecondaryStructureAnnotator.SECONDARY_STRUCTURE_STATES,
            SequenceCerosene.SEQUENCE_CEROSENE_REPRESENTATION,
            SequenceMotifAnnotator.SEQUENCE_MOTIF,
            SiftsParser.UNIPROT_ID,
            PLIPAnnotator.PLIP_INTERACTIONS,
            UniProtAnnotator.UNIPROT_ANNOTATION)
            .map(FeatureProviderRegistry::resolve)
            .collect(Collectors.toList());

    @Autowired
    public ProteinService(ProteinRepository repository) {
        this.repository = repository;
    }

    @PostConstruct
    public void activate() throws IOException {
        allProteins = Stream.of("1m0l").collect(Collectors.toList());

        // clear database
//        repository.deleteAll();
//
//        allProteins.forEach(pdbId -> java.util.concurrent.Executors.newWorkStealingPool().submit(() -> process(pdbId)));
    }

    private void process(String pdbId) {
        try {
            logger.info("fetching information for {}", pdbId);

            // fetch opm and pdb entry
            Protein pdbProtein = ProteinParser.source(pdbId).parse();
            Protein opmProtein = OPMDatabaseQuery.parseAnnotatedProteinById(pdbId);

//            /* the approach where the opm structure is the reference */
//            // superimpose so stuff can be computed on the pdb structure (which is the structure on whose coordinates
//            // PLIP results are available) - do the same for membrane atoms
//            StructureAlignmentResult alignment = new SVDSuperimposer().align(pdbProtein, opmProtein);
//            double[][] rotation = alignment.getRotation();
//            double[] translation = alignment.getTranslation();
//            logger.info("alignment rmsd: {}", alignment.getAlignmentScore());
//            alignment.transform(opmProtein);
//            Membrane membrane = opmProtein.getFeature(Membrane.class, ANVIL.MEMBRANE);
//            List<double[]> transformedMembraneAtoms = membrane.getMembraneAtoms().stream()
//                    .map(coordinates -> LinearAlgebra3D.add(LinearAlgebra3D.multiply(coordinates, rotation), translation))
//                    .collect(Collectors.toList());
//            membrane.setMembraneAtoms(transformedMembraneAtoms);
//
//            // compute features on the
//            opmProtein.setTitle(pdbProtein.getTitle());
//            computeFeatures(opmProtein);
//
//            // persist
//            repository.save(new ExplorerProtein(opmProtein));
//            logger.info("persisted {}", opmProtein.getName());

            /* the approach where the pdb structure is the reference */
            StructureAlignmentResult alignment = new SVDSuperimposer().align(pdbProtein, opmProtein);
            double[][] rotation = alignment.getRotation();
            double[] translation = alignment.getTranslation();
            logger.info("alignment rmsd: {}", alignment.getAlignmentScore());
            Membrane membrane = opmProtein.getFeature(Membrane.class, ANVIL.MEMBRANE);
            List<double[]> transformedMembraneAtoms = membrane.getMembraneAtoms().stream()
                    .map(coordinates -> LinearAlgebra3D.add(LinearAlgebra3D.multiply(coordinates, rotation), translation))
                    .filter(coordinates -> distance(pdbProtein, coordinates))
                    .collect(Collectors.toList());
            membrane.setMembraneAtoms(transformedMembraneAtoms);
            pdbProtein.setFeature(ANVIL.MEMBRANE, membrane);
            pdbProtein.setFeature(OPMDatabaseQuery.HOMOLOGOUS_PROTEINS, opmProtein.getFeatureAsList(String.class, OPMDatabaseQuery.HOMOLOGOUS_PROTEINS));

            // compute features on the pdb protein
            computeFeatures(pdbProtein);

            repository.save(new ExplorerProtein(pdbProtein));
            logger.info("persisted {}", pdbProtein.getName());
        } catch (Exception e) {
            e.printStackTrace();
            logger.error("gathering information for {} failed: {}", pdbId, e.getLocalizedMessage());
        }
    }

    private static boolean distance(Protein protein, final double[] membranePseudoAtom) {
        // too close to protein
        if(protein.atoms()
                .map(Atom::getCoordinates)
                .mapToDouble(coordinates -> LinearAlgebra3D.distanceFast(coordinates, membranePseudoAtom))
                .min()
                .orElse(Double.MAX_VALUE) < 12.0) {
            return false;
        }

        // too far from protein
        if(protein.atoms()
                .map(Atom::getCoordinates)
                .mapToDouble(coordinates -> LinearAlgebra3D.distanceFast(coordinates, membranePseudoAtom))
                .min()
                .orElse(Double.MIN_VALUE) > 90.0) {
            return false;
        }

        //TODO actually compute membrane layers
        return true;
    }

    private void computeFeatures(Protein protein) {
        featureProviders.forEach(featureProvider -> {
            logger.info("processing {} by {}", protein.getName(), featureProvider.getClass().getSimpleName());
            featureProvider.process(protein);
        });
    }

    public List<String> getAllProteins() {
        return allProteins;
    }

    public ExplorerProtein getProtein(String pdbid) {
        List<ExplorerProtein> result = repository.findByTheProteinsName(pdbid.toUpperCase());
        if(result.size() > 0) {
            return result.get(0);
        }

        throw new NoSuchElementException("could not retrieve entry for " + pdbid);
    }

    private Path getResource(String filename) {
        // some a bit hacky way to ensure correct paths on windows (as some / will be added as prefix)
        return Paths.get(Objects.requireNonNull(Thread.currentThread().getContextClassLoader().getResource(filename)).getPath().replaceFirst("^/(.:/)", "$1"));
    }
}
