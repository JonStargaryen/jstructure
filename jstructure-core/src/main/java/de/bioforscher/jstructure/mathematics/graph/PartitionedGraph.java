package de.bioforscher.jstructure.mathematics.graph;

import de.bioforscher.jstructure.mathematics.graph.clustering.Module;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class PartitionedGraph<N> extends Graph<N> {
    private final List<Module<N>> modules;
    private final Module<N> unassignedModule;
    // http://tools.medialab.sciences-po.fr/iwanthue/ - 20 distinct colors
    private static final List<String> PALETTE = Pattern.compile(", ").splitAsStream("[207,61,125], " +
            "[94,198,68], " +
            "[107,61,202], " +
            "[106,156,55], " +
            "[202,74,206], " +
            "[183,171,53], " +
            "[105,111,209], " +
            "[213,136,66], " +
            "[112,46,133], " +
            "[88,181,125], " +
            "[212,73,57], " +
            "[81,169,167], " +
            "[125,47,54], " +
            "[116,148,198], " +
            "[126,85,43], " +
            "[203,122,194], " +
            "[63,95,47], " +
            "[87,64,110], " +
            "[171,165,109], " +
            "[206,134,145]").collect(Collectors.toList());

    public PartitionedGraph(Graph<N> graph, List<Module<N>> modules) {
        super(graph);
        this.modules = modules;
        this.unassignedModule = new Module<>("?", new Graph<N>(new ArrayList<>(), new ArrayList<>()));
    }

    public Module getModuleOf(N node) {
        return modules.stream()
                .filter(module -> module.getNodes().contains(node))
                .findFirst()
                // no module has to associated to a node - if so, return a default value
                .orElse(unassignedModule);
    }

    public List<Module<N>> getModules() {
        return modules;
    }

    public int getNumberOfModules() {
        return modules.size();
    }

    public List<String> getModuleIdentifiers() {
        return modules.stream()
                .map(Module::getIdentifier)
                .collect(Collectors.toList());
    }

    public String getPyMolString() {
        N firstNode = getNodes().get(0);
        //TODO quick-and-dirty implementation with horrible 'design'
        if(!(firstNode instanceof AminoAcid)) {
            throw new UnsupportedOperationException("operation only supported for protein structure graphs");
        }

        AminoAcid aminoAcid = (AminoAcid) firstNode;
        String pdbId = aminoAcid.getParentChain().getChainIdentifier().getProteinIdentifier().getPdbId();
        String chainId = aminoAcid.getParentChain().getChainIdentifier().getChainId();

        return getNodes().stream()
                .map(this::composeColorString)
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.joining(System.lineSeparator(),
                        "fetch " + pdbId + ", async=0" + System.lineSeparator() +
                        "bg_color white" + System.lineSeparator() +
                        // hide non-relevant stuff
                        "hide everything" + System.lineSeparator() +
                        "show cartoon, chain " + chainId + System.lineSeparator() +
                        // decolor everything
                        "color grey80" + System.lineSeparator() +
                        // create custom color palette
                        IntStream.range(0, PALETTE.size())
                                .mapToObj(i -> {
                                    String[] split = PALETTE.get(i)
                                            .replace("[", "")
                                            .replaceAll("]", "")
                                            .split(",");
                                    return "set_color cc" + i + "=[" +
                                            Double.valueOf(split[0]) / 255 + ", " +
                                            Double.valueOf(split[1]) / 255 + ", " +
                                            Double.valueOf(split[2]) / 255 + "]";
                                })
                                .collect(Collectors.joining(System.lineSeparator())) + System.lineSeparator(),
                        ""));
    }

    private Optional<String> composeColorString(N node) {
        Module module = getModuleOf(node);

        // ignore unassigned residues
        if(module.equals(unassignedModule)) {
            return Optional.empty();
        }

        // select and set color
        int moduleNumber = modules.indexOf(module) % PALETTE.size();
        //TODO insertion codes
        int resi = ((AminoAcid) node).getResidueIdentifier().getResidueNumber();

        return Optional.of("color cc" + moduleNumber + ", resi " + resi);
    }
}
