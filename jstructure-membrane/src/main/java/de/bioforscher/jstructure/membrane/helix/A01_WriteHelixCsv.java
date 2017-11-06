package de.bioforscher.jstructure.membrane.helix;

import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.mathematics.Pair;
import de.bioforscher.jstructure.mathematics.SetOperations;
import de.bioforscher.jstructure.membrane.MembraneConstants;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.nio.file.Path;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

/**
 * Evaluates all pairs of TM-helices in the alpha-dataset. Indirectly covers amino acids and contacts.
 */
public class A01_WriteHelixCsv {
    private static final Logger logger = LoggerFactory.getLogger(A01_WriteHelixCsv.class);
    private static final Path directory = MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY;

    public static void main(String[] args) {
        // general stuff
        String header = "pdbId;chainId;tm1;tm2;start1;end1;start2;end2;" + System.lineSeparator();

        String output = MembraneConstants.lines(directory.resolve("ids.list"))
                .map(A01_WriteHelixCsv::handleLine)
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.joining(System.lineSeparator()));

        MembraneConstants.write(directory.resolve("results").resolve("helices.csv"),
                header + output);
    }

    private static Optional<String> handleLine(String line) {
        try {
            Optional<MembraneConstants.WrappedChain> wrappedChainOptional = MembraneConstants.handleLine(line, directory);

            if(!wrappedChainOptional.isPresent()) {
                return Optional.empty();
            }

            MembraneConstants.WrappedChain wrappedChain = wrappedChainOptional.get();
            Chain chain = wrappedChain.getChain();
            String pdbId = wrappedChain.getPdbId();
            String chainId = wrappedChain.getChainId();
            boolean evolScoresSane = wrappedChain.isEvolScoresSane();
            boolean kinkFinderDataSane = wrappedChain.isKinkFinderDataSane();
            List<List<AminoAcid>> transmembraneHelices = wrappedChain.getTransmembraneHelices();


            return Optional.of(SetOperations.uniquePairsOf(transmembraneHelices)
                    .map(pair -> handleHelixPair(pair,
                            transmembraneHelices,
                            pdbId,
                            chainId,
                            chain,
                            evolScoresSane,
                            kinkFinderDataSane))
                    .filter(Optional::isPresent)
                    .map(Optional::get)
                    .collect(Collectors.joining(System.lineSeparator())));
        } catch (Exception e) {
            e.printStackTrace();
            return Optional.empty();
        }
    }

    private static Optional<String> handleHelixPair(Pair<List<AminoAcid>, List<AminoAcid>> pair,
                                                    List<List<AminoAcid>> transmembraneHelices,
                                                    String pdbId,
                                                    String chainId,
                                                    Chain chain,
                                                    boolean evolScoresSane,
                                                    boolean kinkFinderDataSane) {
        try {
            int index1 = transmembraneHelices.indexOf(pair.getLeft()) + 1;
            int index2 = transmembraneHelices.indexOf(pair.getRight()) + 1;

            PLIPInteractionContainer contacts = new PLIPInteractionContainer(null,
                    chain.getFeature(PLIPInteractionContainer.class)
                            .getInteractions()
                            .stream()
                            .filter(plipInteraction -> (pair.getLeft().contains(plipInteraction.getPartner1()) && pair.getRight().contains(plipInteraction.getPartner2())) ||
                                    (pair.getLeft().contains(plipInteraction.getPartner2()) && pair.getRight().contains(plipInteraction.getPartner1())))
                            .collect(Collectors.toList()));

            logger.info("{} contact between helix {} and {}",
                    contacts.getInteractions().size(),
                    index1,
                    index2);

            return Optional.of(pdbId + ";" +
                    chainId + ";" +
                    index1 + ";" +
                    index2 + ";" +
                    pair.getLeft().get(0).getResidueIdentifier() + ";" +
                    pair.getLeft().get(pair.getLeft().size() - 1).getResidueIdentifier() + ";" +
                    pair.getRight().get(0).getResidueIdentifier() + ";" +
                    pair.getRight().get(pair.getRight().size() - 1).getResidueIdentifier() + ";" +
                    contacts.getHalogenBonds().size() + ";" +
                    contacts.getHydrogenBonds().size() + ";" +
                    contacts.getHydrophobicInteractions().size() + ";" +
                    contacts.getMetalComplexes().size() + ";" +
                    contacts.getPiCationInteractions().size() + ";" +
                    contacts.getPiStackings().size() + ";" +
                    contacts.getSaltBridges().size() + ";" +
                    contacts.getWaterBridges().size() + ";" +
                    contacts.getInteractions().size()
            );
        } catch (Exception e) {
            logger.warn("something something", e);
            return Optional.empty();
        }
    }
}
