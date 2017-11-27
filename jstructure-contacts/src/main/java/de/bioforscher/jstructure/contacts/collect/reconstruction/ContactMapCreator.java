package de.bioforscher.jstructure.contacts.collect.reconstruction;

import de.bioforscher.jstructure.contacts.ContactsConstants;
import de.bioforscher.jstructure.feature.interactions.*;
import de.bioforscher.jstructure.mathematics.Pair;
import de.bioforscher.jstructure.mathematics.SetOperations;
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

public class ContactMapCreator {
    private static final Logger logger = LoggerFactory.getLogger(ContactMapCreator.class);

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

            logger.info("sampling contact list with {} contacts - keeping {}%, i.e. {} contacts",
                    size,
                    fraction,
                    keeping);

            return keeping;
        }

        public List<Edge<AminoAcid>> sampleContactListRandomly(List<Edge<AminoAcid>> edges) {
            int keeping = reportKeeps(edges);

            Collections.shuffle(edges);

            return edges.subList(0, keeping);
        }
    }

    public ContactMapCreator(Chain chain, Path contactMapDirectory, Path start2foldPath) {
        List<Group> earlyFoldingResidues = ContactsConstants.lines(start2foldPath)
                .filter(line -> !line.startsWith("#"))
                .filter(line -> line.endsWith("EARLY"))
                .map(line -> line.split(";")[0])
                .map(Integer::valueOf)
                .map(residueNumber -> chain.select().residueNumber(residueNumber).asOptionalAminoAcid())
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.toList());
        List<Edge<AminoAcid>> averagedPlipEdges = chain.getFeature(PLIPInteractionContainer.class)
                .getInteractions()
                .stream()
                // ignore covalently bound groups
                .filter(plipInteraction -> areNonCovalentGroups(plipInteraction.getPartner1(), plipInteraction.getPartner2()))
                .filter(plipInteraction -> isEarlyFoldingInteraction(plipInteraction, earlyFoldingResidues))
                .map(plipInteraction -> new Edge<>(((AminoAcid) plipInteraction.getPartner1()), ((AminoAcid) plipInteraction.getPartner2()), determinePlipInteractionDistanceByAverage(plipInteraction)))
                .collect(Collectors.toList());
        List<Edge<AminoAcid>> maximumPlipEdges = chain.getFeature(PLIPInteractionContainer.class)
                .getInteractions()
                .stream()
                // ignore covalently bound groups
                .filter(plipInteraction -> areNonCovalentGroups(plipInteraction.getPartner1(), plipInteraction.getPartner2()))
                .filter(plipInteraction -> isEarlyFoldingInteraction(plipInteraction, earlyFoldingResidues))
                .map(plipInteraction -> new Edge<>(((AminoAcid) plipInteraction.getPartner1()), ((AminoAcid) plipInteraction.getPartner2()), determinePlipInteractionDistanceByMaximum(plipInteraction)))
                .collect(Collectors.toList());

        int plipEdgeCount = averagedPlipEdges.size();
        List<Edge<AminoAcid>> naiveEdges = determineNaiveInteractions(chain.aminoAcids().collect(Collectors.toList()));

        String id = chain.getChainIdentifier().getFullName();
        String sequence = chain.getAminoAcidSequence();

        ContactsConstants.write(contactMapDirectory.resolve(id + ".fasta"),
                ">" + id + System.lineSeparator() + sequence);

        ContactsConstants.write(contactMapDirectory.resolve(id + "-early-avg.rr"),
                composeRRString(averagedPlipEdges, sequence));
        ContactsConstants.write(contactMapDirectory.resolve(id + "-early-max.rr"),
                composeRRString(averagedPlipEdges, sequence));
        // create 10 sampled maps with an equal number of contacts
        for(int i = 1; i < 11; i++) {
            Collections.shuffle(naiveEdges);

            List<Edge<AminoAcid>> sampledNaiveEdges = naiveEdges.subList(0, plipEdgeCount);
            ContactsConstants.write(contactMapDirectory.resolve(id + "-naive-" + i + ".rr"),
                    composeRRString(sampledNaiveEdges, sequence));
        }
    }

    public ContactMapCreator(Chain chain, Path contactMapDirectory) {
        List<Edge<AminoAcid>> plipEdges = chain.getFeature(PLIPInteractionContainer.class)
                .getInteractions()
                .stream()
                // ignore covalently bound groups
                .filter(plipInteraction -> areNonCovalentGroups(plipInteraction.getPartner1(), plipInteraction.getPartner2()))
                .map(plipInteraction -> new Edge<>(((AminoAcid) plipInteraction.getPartner1()), ((AminoAcid) plipInteraction.getPartner2()), determinePlipInteractionDistanceByAverage(plipInteraction)))
                .collect(Collectors.toList());
        List<Edge<AminoAcid>> naiveEdges = determineNaiveInteractions(chain.aminoAcids().collect(Collectors.toList()));

        String id = chain.getChainIdentifier().getFullName();
        String sequence = chain.getAminoAcidSequence();

        ContactsConstants.write(contactMapDirectory.resolve(id + ".fasta"),
                ">" + id + System.lineSeparator() + sequence);

        ContactsConstants.write(contactMapDirectory.resolve("p100").resolve(id + "-plip.rr"),
                composeRRString(plipEdges, sequence));
        ContactsConstants.write(contactMapDirectory.resolve("p100").resolve(id + "-naive.rr"),
                composeRRString(naiveEdges, sequence));

        for(Sampling sampling : Sampling.values()) {
            for(int i = 0; i < SAMPLING_SIZE; i++) {
                List<Edge<AminoAcid>> plipEdgesSampled = sampling.sampleContactListRandomly(plipEdges);
                List<Edge<AminoAcid>> naiveEdgesSampled = sampling.sampleContactListRandomly(naiveEdges);

                ContactsConstants.write(contactMapDirectory.resolve(sampling.name()).resolve(id + "-plip-" + (i + 1) + ".rr"),
                        composeRRString(plipEdgesSampled, sequence));
                ContactsConstants.write(contactMapDirectory.resolve(sampling.name()).resolve(id + "-naive-" + (i + 1) + ".rr"),
                        composeRRString(naiveEdgesSampled, sequence));
            }
        }
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
    private static final double ALPHA_THRESHOLD = 9.0;
    private static final double BETA_THRESHOLD = 8.0;

    private List<Edge<AminoAcid>> determineNaiveInteractions(List<AminoAcid> aminoAcids) {
        List<Edge<AminoAcid>> naiveEdges = new ArrayList<>();
        List<Pair<AminoAcid, AminoAcid>> pairs = SetOperations.uniquePairsOf(aminoAcids).collect(Collectors.toList());

        for(Pair<AminoAcid, AminoAcid> pair : pairs) {
            AminoAcid aminoAcid1 = pair.getLeft();
            AminoAcid aminoAcid2 = pair.getRight();

            if(!areNonCovalentGroups(aminoAcid1, aminoAcid2)) {
                continue;
            }

//            if(aminoAcid1.getCa().calculate().distance(aminoAcid2.getCa()) < ALPHA_THRESHOLD) {
//                naiveEdges.add(new Edge<>(aminoAcid1, aminoAcid2, ALPHA_THRESHOLD));
//            }

            if(ContactsConstants.getBetaCarbon(aminoAcid1).calculate().distance(ContactsConstants.getBetaCarbon(aminoAcid2)) < BETA_THRESHOLD) {
                naiveEdges.add(new Edge<>(aminoAcid1, aminoAcid2, BETA_THRESHOLD));
            }
        }

        return naiveEdges;
    }

    private double determinePlipInteractionDistanceByMaximum(PLIPInteraction plipInteraction) {
        // max values
        if(plipInteraction instanceof SaltBridge) {
            return 10.8307;
        } else if(plipInteraction instanceof WaterBridge) {
            return 12.1106;
        } else if(plipInteraction instanceof PiCationInteraction) {
            return 11.9862;
        } else if(plipInteraction instanceof HydrophobicInteraction) {
            return 5.0976;
        } else if(plipInteraction instanceof PiStacking) {
            return 7.393;
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

    private boolean isEarlyFoldingInteraction(PLIPInteraction plipInteraction, List<Group> earlyFoldingResidues) {
        return earlyFoldingResidues.contains(plipInteraction.getPartner1()) || earlyFoldingResidues.contains(plipInteraction.getPartner2());
    }


    private boolean areNonCovalentGroups(Group group1, Group group2) {
        return Math.abs(group1.getResidueIndex() - group2.getResidueIndex()) > 1;
    }
}
