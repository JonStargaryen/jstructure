package de.bioforscher.jstructure.membrane.modularity.pdbtm;

import de.bioforscher.jstructure.membrane.MembraneConstants;
import de.bioforscher.jstructure.membrane.modularity.GraphFactory;
import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.SetOperations;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.nio.file.Path;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

/**
 * Execute the NetCarlo-algorithm on a dataset.
 * TODO can multiple links between nodes emulate a connections strength? i.e. hbonds > hydrophobic interactions
 */
public class NetworkExtractor {
    private static final Logger logger = LoggerFactory.getLogger(NetworkExtractor.class);
    private static final double STRICT_NAIVE_DISTANCE_CUTOFF = 0.6;
    private static final double GENEROUS_NAIVE_DISTANCE_CUTOFF = 1.0;
    private final Path pdbPath;
    private final Path plipPath;
    private final Path networkPath;

    public NetworkExtractor(Path datasetPath) {
        Path listPath = datasetPath.resolve("ids.list");
        this.pdbPath = datasetPath.resolve("pdb");
        this.plipPath = datasetPath.resolve("plip");
        this.networkPath = datasetPath.resolve("network");

        MembraneConstants.lines(listPath)
                .forEach(this::handleLine);
    }

    private void handleLine(String id) {
        logger.info("processing {}", id);

        String pdbId = id.split("_")[0];
        String chainId = id.split("_")[1];

        Chain chain = StructureParser.source(pdbPath.resolve(pdbId + ".pdb"))
                .minimalParsing(true)
                .parse()
                .select()
                .chainId(chainId)
                .asChain();

        // extract network files
        Pair<String, String> naiveNetworks = getNaiveNetwork(chain);
//        MembraneConstants.write(networkPath.resolve(id + "_strict_nc.dat"), naiveNetworks.getLeft());
//        MembraneConstants.write(networkPath.resolve(id + "_generous_nc.dat"), naiveNetworks.getRight());
        MembraneConstants.write(networkPath.resolve(id + "_strict.dat"), naiveNetworks.getLeft());
        MembraneConstants.write(networkPath.resolve(id + "_generous.dat"), naiveNetworks.getRight());

        String document = MembraneConstants.lines(plipPath.resolve(id + ".plip")).collect(Collectors.joining(System.lineSeparator()));
        MembraneConstants.write(networkPath.resolve(id + "_plip.dat"), getPlipNetwork(chain, document));
//        MembraneConstants.write(networkPath.resolve(id + "_plip_nc.dat"), getPlipNetwork(chain, document));
    }

    private String getPlipNetwork(Chain chain, String document) {
        return GraphFactory.createGraphFromPlipDocument(chain, document).toNetCartoString();
    }

    /**
     * Defines contacts of amino acids as pairs who share a pair of atoms less than a defined threshold apart.
     * @param chain the structure to process
     * @return a string in NetCarlo's input format of observed interactions/contacts
     */
    private Pair<String, String> getNaiveNetwork(Chain chain) {
        List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());
        StringJoiner strictContacts = new StringJoiner(System.lineSeparator());
        StringJoiner generousContacts = new StringJoiner(System.lineSeparator());

        // russian way of doing both thresholds in the same run as calculation is somewhat taxing
        for(int i = 0; i < aminoAcids.size() - 1; i++) {
            AminoAcid aminoAcid1 = aminoAcids.get(i);
            List<Atom> atoms1 = aminoAcid1.getAtoms();
            for(int j = i + 1; j < aminoAcids.size(); j++) {
                AminoAcid aminoAcid2 = aminoAcids.get(j);
                if(j == i + 1) {
                    // comment next two lines to ignore consecutive amino acids
                    generousContacts.add(aminoAcid1.getResidueIdentifier() + " " + aminoAcid2.getResidueIdentifier());
                    strictContacts.add(aminoAcid1.getResidueIdentifier() + " " + aminoAcid2.getResidueIdentifier());
                    continue;
                }

                List<Atom> atoms2 = aminoAcid2.getAtoms();
                if(SetOperations.cartesianProductOf(atoms1, atoms2).anyMatch(this::isGenerousNaiveContact)) {
                    generousContacts.add(aminoAcid1.getResidueIdentifier() + " " + aminoAcid2.getResidueIdentifier());
                } else {
                    // if generous threshold is not met, the strict one can neither
                    continue;
                }
                if(SetOperations.cartesianProductOf(atoms1, atoms2).anyMatch(this::isStrictNaiveContact)) {
                    strictContacts.add(aminoAcid1.getResidueIdentifier() + " " + aminoAcid2.getResidueIdentifier());
                }
            }
        }

        return new Pair<>(strictContacts.toString(), generousContacts.toString());
    }

    private boolean isStrictNaiveContact(Pair<Atom, Atom> pair) {
        return isNaiveContact(pair, STRICT_NAIVE_DISTANCE_CUTOFF);
    }

    private boolean isGenerousNaiveContact(Pair<Atom, Atom> pair) {
        return isNaiveContact(pair, GENEROUS_NAIVE_DISTANCE_CUTOFF);
    }

    private boolean isNaiveContact(Pair<Atom, Atom> pair, double distanceCutoff) {
        Atom a1 = pair.getLeft();
        Atom a2 = pair.getRight();
        double distance = a1.calculate().distance(a2);
        return distance <= a1.getElement().getVDWRadius() + a2.getElement().getVDWRadius() + distanceCutoff;
    }
}
