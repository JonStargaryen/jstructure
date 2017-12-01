package de.bioforscher.jstructure.contacts.collect.reconstruction;

import de.bioforscher.jstructure.contacts.ContactsConstants;
import de.bioforscher.jstructure.feature.interactions.*;
import de.bioforscher.jstructure.mathematics.graph.Edge;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

public class ReconstructionContactMapCreator {
    private static final Logger logger = LoggerFactory.getLogger(ReconstructionContactMapCreator.class);

    public enum Sampling {
        p10(10),
        p30(30),
        p50(50),
        p70(70),
        p90(90);

        private final int fraction;

        Sampling(int fraction) {
            this.fraction = fraction;
        }

        public int getFraction() {
            return fraction;
        }

        private int reportKeeps(List<Edge<AminoAcid>> edges) {
            int size = edges.size();
            int keeping = (int) (size * (((double) fraction) / 100.0));

            // ensure at least 1 contact
            if(keeping == 0) {
                keeping = 1;
            }

//            logger.info("sampling contact list with {} contacts - keeping {}%, i.e. {} contacts",
//                    size,
//                    fraction,
//                    keeping);

            return keeping;
        }

        public List<Edge<AminoAcid>> sampleContactListRandomly(List<Edge<AminoAcid>> edges) {
            int keeping = reportKeeps(edges);

            Collections.shuffle(edges);

            return edges.subList(0, keeping);
        }
    }

    public ReconstructionContactMapCreator(Chain chain, Path contactMapDirectory, Path start2foldPath) {
        List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());
        List<Group> earlyFoldingResidues = ContactsConstants.lines(start2foldPath)
                .filter(line -> !line.startsWith("#"))
                .filter(line -> line.endsWith("EARLY"))
                .map(line -> line.split(";")[0])
                .map(Integer::valueOf)
                .map(residueNumber -> chain.select().residueNumber(residueNumber).asOptionalAminoAcid())
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.toList());

        // handle PLIP interaction scheme
        List<Edge<AminoAcid>> plipEdges = chain.getFeature(PLIPInteractionContainer.class)
                .getInteractions()
                .stream()
                // ignore covalently bound groups
                .filter(plipInteraction -> ContactsConstants.areNonCovalentGroups(plipInteraction.getPartner1(), plipInteraction.getPartner2()))
                .map(plipInteraction -> new Edge<>(((AminoAcid) plipInteraction.getPartner1()), ((AminoAcid) plipInteraction.getPartner2()), determinePlipInteractionDistanceByMaximum(plipInteraction)))
                .collect(Collectors.toList());
        List<Edge<AminoAcid>> earlyFoldingPlipEdges = plipEdges.stream()
                .filter(edge -> isEarlyFoldingEdge(earlyFoldingResidues, edge))
                .collect(Collectors.toList());
        List<Edge<AminoAcid>> lateFoldingPlipEdges = plipEdges.stream()
                .filter(edge -> !isEarlyFoldingEdge(earlyFoldingResidues, edge))
                .collect(Collectors.toList());

        int earlyFoldingPlipEdgeCount = earlyFoldingPlipEdges.size();
        int lateFoldingPlipEdgeCount = lateFoldingPlipEdges.size();

        String id = chain.getChainIdentifier().getFullName();
        String sequence = chain.getAminoAcidSequence();

        ContactsConstants.write(contactMapDirectory.resolve(id + ".fasta"),
                ">" + id + System.lineSeparator() + sequence);

        ContactsConstants.write(contactMapDirectory.resolve(id + "-early-plip.rr"),
                composeRRString(earlyFoldingPlipEdges, sequence));
        ContactsConstants.write(contactMapDirectory.resolve(id + "-late-plip.rr"),
                composeRRString(lateFoldingPlipEdges, sequence));
        // create 10 sampled maps with an equal number of contacts
        for(int i = 1; i < 11; i++) {
            Collections.shuffle(plipEdges);

            List<Edge<AminoAcid>> sampledEarlyPlipEdges = plipEdges.subList(0, earlyFoldingPlipEdgeCount);
            ContactsConstants.write(contactMapDirectory.resolve(id + "-sampled-plip-" + i + ".rr"),
                    composeRRString(sampledEarlyPlipEdges, sequence));

            List<Edge<AminoAcid>> selectedEarlyPlipEdges = selectEdges(plipEdges, aminoAcids, earlyFoldingPlipEdgeCount);
            ContactsConstants.write(contactMapDirectory.resolve(id + "-selected-plip-" + i + ".rr"),
                    composeRRString(selectedEarlyPlipEdges, sequence));

            List<Edge<AminoAcid>> sampledLatePlipEdges = plipEdges.subList(0, lateFoldingPlipEdgeCount);
            ContactsConstants.write(contactMapDirectory.resolve(id + "-late_sampled-plip-" + i + ".rr"),
                    composeRRString(sampledLatePlipEdges, sequence));

            List<Edge<AminoAcid>> selectedLatePlipEdges = selectEdges(plipEdges, aminoAcids, lateFoldingPlipEdgeCount);
            ContactsConstants.write(contactMapDirectory.resolve(id + "-late_selected-plip-" + i + ".rr"),
                    composeRRString(selectedLatePlipEdges, sequence));
        }

        // handle conventional interaction scheme
        List<Edge<AminoAcid>> conventionalEdges = ContactsConstants.determineNaiveInteractions(chain.aminoAcids().collect(Collectors.toList()));
        List<Edge<AminoAcid>> earlyFoldingConventionalEdges = conventionalEdges.stream()
                .filter(edge -> isEarlyFoldingEdge(earlyFoldingResidues, edge))
                .collect(Collectors.toList());
        List<Edge<AminoAcid>> lateFoldingConventionalEdges = conventionalEdges.stream()
                .filter(edge -> !isEarlyFoldingEdge(earlyFoldingResidues, edge))
                .collect(Collectors.toList());

        int earlyFoldingConventionalEdgeCount = earlyFoldingConventionalEdges.size();
        int lateFoldingConventionalEdgeCount = lateFoldingConventionalEdges.size();

        ContactsConstants.write(contactMapDirectory.resolve(id + "-early-conventional.rr"),
                composeRRString(earlyFoldingConventionalEdges, sequence));
        ContactsConstants.write(contactMapDirectory.resolve(id + "-late-conventional.rr"),
                composeRRString(lateFoldingConventionalEdges, sequence));
        // create 10 sampled maps with an equal number of contacts
        for(int i = 1; i < 11; i++) {
            Collections.shuffle(conventionalEdges);

            List<Edge<AminoAcid>> sampledEarlyConventionalEdges = conventionalEdges.subList(0, earlyFoldingConventionalEdgeCount);
            ContactsConstants.write(contactMapDirectory.resolve(id + "-sampled-conventional-" + i + ".rr"),
                    composeRRString(sampledEarlyConventionalEdges, sequence));

            List<Edge<AminoAcid>> selectedEarlyConventionalEdges = selectEdges(conventionalEdges, aminoAcids, earlyFoldingConventionalEdgeCount);
            ContactsConstants.write(contactMapDirectory.resolve(id + "-selected-conventional-" + i + ".rr"),
                    composeRRString(selectedEarlyConventionalEdges, sequence));

            List<Edge<AminoAcid>> sampledLateConventionalEdges = conventionalEdges.subList(0, lateFoldingConventionalEdgeCount);
            ContactsConstants.write(contactMapDirectory.resolve(id + "-late_sampled-conventional-" + i + ".rr"),
                    composeRRString(sampledLateConventionalEdges, sequence));

            List<Edge<AminoAcid>> selectedLateConventionalEdges = selectEdges(conventionalEdges, aminoAcids, lateFoldingConventionalEdgeCount);
            ContactsConstants.write(contactMapDirectory.resolve(id + "-late_selected-conventional-" + i + ".rr"),
                    composeRRString(selectedLateConventionalEdges, sequence));
        }
    }

    private List<Edge<AminoAcid>> selectEdges(List<Edge<AminoAcid>> edges, List<AminoAcid> aminoAcids, int count) {
        Collections.shuffle(aminoAcids);
        List<Edge<AminoAcid>> selectedEdges = new ArrayList<>();

        //TODO here as well as with the other selections an edge may get selected twice
        for(AminoAcid aminoAcid : aminoAcids) {
            edges.stream()
                    .filter(edge ->  edge.contains(aminoAcid))
                    .forEach(selectedEdges::add);

            // select at least as many edges as there are early folding ones
            if(selectedEdges.size() > count) {
                break;
            }
        }

        return selectedEdges;
    }

    public ReconstructionContactMapCreator(Chain chain, Path contactMapDirectory) {
        List<Edge<AminoAcid>> plipEdges = chain.getFeature(PLIPInteractionContainer.class)
                .getInteractions()
                .stream()
                // ignore covalently bound groups
                .filter(plipInteraction -> ContactsConstants.areNonCovalentGroups(plipInteraction.getPartner1(), plipInteraction.getPartner2()))
                .map(plipInteraction -> new Edge<>(((AminoAcid) plipInteraction.getPartner1()), ((AminoAcid) plipInteraction.getPartner2()), determinePlipInteractionDistanceByMaximum(plipInteraction)))
                .collect(Collectors.toList());
        List<Edge<AminoAcid>> naiveEdges = ContactsConstants.determineNaiveInteractions(chain.aminoAcids().collect(Collectors.toList()));
//        List<Edge<AminoAcid>> enrichedEdges = enrichNaiveInteractions(naiveEdges, plipEdges);

        String id = chain.getChainIdentifier().getFullName();
        String sequence = chain.getAminoAcidSequence();

        ContactsConstants.write(contactMapDirectory.resolve(id + ".fasta"),
                ">" + id + System.lineSeparator() + sequence);

        ContactsConstants.write(contactMapDirectory.resolve("p100").resolve(id + "-plip.rr"),
                composeRRString(plipEdges, sequence));
        ContactsConstants.write(contactMapDirectory.resolve("p100").resolve(id + "-naive.rr"),
                composeRRString(naiveEdges, sequence));
//        ContactsConstants.write(contactMapDirectory.resolve("p100").resolve(id + "-enrich.rr"),
//                composeRRString(enrichedEdges, sequence));

        for(Sampling sampling : Sampling.values()) {
            for(int i = 0; i < SAMPLING_SIZE; i++) {
                List<Edge<AminoAcid>> plipEdgesSampled = sampling.sampleContactListRandomly(plipEdges);
                List<Edge<AminoAcid>> naiveEdgesSampled = sampling.sampleContactListRandomly(naiveEdges);
//                List<Edge<AminoAcid>> enrichedEdgesSampled = sampling.sampleContactListRandomly(enrichedEdges);
//
                ContactsConstants.write(contactMapDirectory.resolve(sampling.name()).resolve(id + "-plip-" + (i + 1) + ".rr"),
                        composeRRString(plipEdgesSampled, sequence));
                ContactsConstants.write(contactMapDirectory.resolve(sampling.name()).resolve(id + "-naive-" + (i + 1) + ".rr"),
                        composeRRString(naiveEdgesSampled, sequence));
//                ContactsConstants.write(contactMapDirectory.resolve(sampling.name()).resolve(id + "-enrich-" + (i + 1) + ".rr"),
//                        composeRRString(enrichedEdgesSampled, sequence));
            }
        }
    }

    /**
     * Merge both contact definitions and employ more strict upper bounds for shared interactions.
     * performance:
     * - rejected: occurs for about 5 constraints in the whole dataset
     * @param naiveEdges
     * @param plipEdges
     * @return the set of naive edges but with more strict distance constraints
     */
    private List<Edge<AminoAcid>> enrichNaiveInteractions(List<Edge<AminoAcid>> naiveEdges, List<Edge<AminoAcid>> plipEdges) {
        return naiveEdges.stream()
                .map(naiveEdge -> {
                    double weight = plipEdges.stream()
                            .filter(plipEdge -> plipEdge.contains(naiveEdge.getLeft()) && plipEdge.contains(naiveEdge.getRight()))
                            .mapToDouble(Edge::getWeight)
                            .filter(distance -> distance < BETA_THRESHOLD)
                            // select strongest constraint if there are multiple
                            .min()
                            // if none is present: keep weight of initial
                            .orElse(BETA_THRESHOLD);
                    if(weight < BETA_THRESHOLD) {
                        logger.info("refining constraint between {} and {} from {} to {}",
                                naiveEdge.getLeft(),
                                naiveEdge.getRight(),
                                BETA_THRESHOLD,
                                weight);
                    }
                    return new Edge<>(naiveEdge.getLeft(), naiveEdge.getRight(), weight);
                })
                .collect(Collectors.toList());
    }

    private String composeRRString(List<Edge<AminoAcid>> edges, String sequence) {
        // i  j  d1  d2  p
        // 24 33 0 8 12.279515
        // i and j: residue numbers
        return edges.stream()
                .map(edge -> (edge.getLeft().getResidueIndex() + 1) + " " +
                        (edge.getRight().getResidueIndex() + 1) + " " +
                        "0 " +
                        edge.getWeight() + " " +
                        1)
                .collect(Collectors.joining(System.lineSeparator(),
                        sequence + System.lineSeparator(),
                        ""));
    }

    private static final int SAMPLING_SIZE = 10;
    private static final double BETA_THRESHOLD = 8.0;

    private double determinePlipInteractionDistanceByMaximum(PLIPInteraction plipInteraction) {
        // max values - TODO move to nr-PDB dataset?
        if(plipInteraction instanceof SaltBridge) {
            return 10.8307;
        } else if(plipInteraction instanceof WaterBridge) {
            return 12.1106;
        } else if(plipInteraction instanceof PiCationInteraction) {
            return 7.393;
        } else if(plipInteraction instanceof HydrophobicInteraction) {
            return 11.9862;
        } else if(plipInteraction instanceof PiStacking) {
            return 5.0976;
        } else if(plipInteraction instanceof HydrogenBond) {
            return 11.5943;
        } else {
            return 12.1106;
        }
    }

    private double determinePlipInteractionDistanceByAverage(PLIPInteraction plipInteraction) {
        // average distance values determine on this particular data - TODO move to nr-PDB dataset?
        if(plipInteraction instanceof SaltBridge) {
            return 6.8865;
        } else if(plipInteraction instanceof WaterBridge) {
            return 6.5379;
        } else if(plipInteraction instanceof PiCationInteraction) {
            return 5.8719;
        } else if(plipInteraction instanceof HydrophobicInteraction) {
            return 5.9347;
        } else if(plipInteraction instanceof PiStacking) {
            return 4.7513;
        } else if(plipInteraction instanceof HydrogenBond) {
            return 5.8208;
        } else {
            return 5.9308;
        }
    }

    private boolean isEarlyFoldingEdge(List<Group> earlyFoldingResdiues, Edge<AminoAcid> edge) {
        return earlyFoldingResdiues.contains(edge.getLeft()) || earlyFoldingResdiues.contains(edge.getRight());
    }

    private boolean isEarlyFoldingInteraction(PLIPInteraction plipInteraction, List<Group> earlyFoldingResidues) {
        return earlyFoldingResidues.contains(plipInteraction.getPartner1()) || earlyFoldingResidues.contains(plipInteraction.getPartner2());
    }
}
