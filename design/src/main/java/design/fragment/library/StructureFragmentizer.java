package design.fragment.library;

import de.bioforscher.jstructure.alignment.consensus.FragmentClusteringComposer;
import de.bioforscher.jstructure.alignment.consensus.StructureCluster;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.selection.Selection;

import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * A delegating wrapper for {@link FragmentClusteringComposer}
 * Created by S on 29.12.2016.
 */
public class StructureFragmentizer {
    /**
     * The fragment size into which protein structures are split. 5 seems reasonable as it is in agreement with the most
     * popular Gerstein motifs such as GG4 which describes 1 turn of a TM-helix.
     */
    private static final int DEFAULT_FRAGMENT_SIZE = 5;

    /**
     * The cutoff below 2 fragments are considered to be similar.
     */
    private static final double DEFAULT_RMSD_CUTOFF = 1.0;

    private final FragmentClusteringComposer fragmentClusteringComposer;

    public StructureFragmentizer() {
        this.fragmentClusteringComposer = new FragmentClusteringComposer(DEFAULT_RMSD_CUTOFF, false);
    }

    public List<StructureCluster> fragmentize(List<Protein> proteins) {
        List<AtomContainer> fragments = proteins.parallelStream()
                .map(this::fragmentize)
                .flatMap(Collection::stream)
                .collect(Collectors.toList());

        System.out.println("saw " + fragments.size() + " fragments");

        fragmentClusteringComposer.composeClusterRepresentation(fragments);

        return fragmentClusteringComposer.getClusters();
    }

    private List<AtomContainer> fragmentize(Protein protein) {
        return protein.getChains().stream()
                .flatMap(this::fragmentize)
                .collect(Collectors.toList());
    }

    private Stream<AtomContainer> fragmentize(Chain chain) {
        List<Group> aminoAcids = Selection.on(chain)
                .aminoAcids()
                .asGroupContainer()
                .getGroups();
        return IntStream.range(0, aminoAcids.size() - DEFAULT_FRAGMENT_SIZE + 1)
                .mapToObj(i -> {
                    Chain container = new Chain(aminoAcids.subList(i, i + DEFAULT_FRAGMENT_SIZE));
                    // need to set identifier of fragments
                    //TODO move to Selection-API?
                    container.setIdentifier(chain.getParentProtein().getIdentifier() + "-" + chain.getChainId() + "-" + i + "-" + (i + DEFAULT_FRAGMENT_SIZE));
                    return container;
                });
    }
}