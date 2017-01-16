package design.fragment.library;

import de.bioforscher.jstructure.model.Combinatorics;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import design.DesignConstants;
import design.ProteinSource;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Move fragments into clusters by the naive approach, but this time not evaluated by the RMSD, now using torsion angles
 * instead.
 * Created by bittrich on 1/9/17.
 */
public class S02_BuildLibraryByTorsionHashing {
    private static final String basePath = DesignConstants.NAIVE_FRAGMENT_CLUSTERS_BY_TORSION_HASHING_DIR;
    private final List<GroupContainer> fragments;
    private final Map<int[], List<AtomContainer>> clusters;

    public void main(String[] args) throws IOException {
        run();
    }

    private static void run() {
        S02_BuildLibraryByTorsionHashing s02_buildLibraryByTorsionHashing = new S02_BuildLibraryByTorsionHashing();
        DesignConstants.makeDirectoryIfAbsent(Paths.get(basePath));
        DesignConstants.TOPOLOGIES.forEach(s02_buildLibraryByTorsionHashing::handleTopology);
    }

    private S02_BuildLibraryByTorsionHashing() {
        this.clusters = new HashMap<>();
        // load proteins
        List<Protein> proteins = ProteinSource.loadProteins(false, false, true);
        this.fragments = new ArrayList<>();

        for(Protein protein : proteins) {
            for(Chain chain : protein.getChains()) {
                List<Group> groups = Selection.on(chain).aminoAcids().asFilteredGroups().collect(Collectors.toList();
                if(groups.size() < 5) {
                    continue;
                }

                Combinatorics.fragmentsOf(groups, 5).forEach(fragments::add);
            }
        }
    }

    private void handleTopology(String topology) {
        String topologyPath = basePath + topology + "/";
        DesignConstants.makeDirectoryIfAbsent(Paths.get(topology));

        // split everything into fragments

    }
}
