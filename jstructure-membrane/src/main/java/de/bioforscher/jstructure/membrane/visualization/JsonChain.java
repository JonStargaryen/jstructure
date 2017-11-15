package de.bioforscher.jstructure.membrane.visualization;

import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceArea;
import de.bioforscher.jstructure.feature.energyprofile.EnergyProfile;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.feature.topology.Topology;
import de.bioforscher.jstructure.mathematics.graph.PartitionedGraph;
import de.bioforscher.jstructure.mathematics.graph.partitioning.Module;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.apache.commons.math3.stat.StatUtils;

import java.util.*;
import java.util.stream.Collectors;

public class JsonChain {
    private final String sequence;
    private final List<Integer> early;
    private final List<JsonSecondaryStructure> sse;
    private final List<Integer> interactions;
    private List<Double> energy;
    private final List<Double> rasa;
    private final List<Double> dynamine;
    private final List<Double> efoldmine;
    private final List<Double> z;
    private final List<Double> p;
    private final List<JsonPartitioning> partitionings;
    private final List<JsonRange> topology;

    public JsonChain(Chain chain, Map<String, PartitionedGraph<AminoAcid>> inSilicoData) {
        this.sequence = chain.getAminoAcidSequence();
        this.sse = composeSecondaryStructures(chain);

        this.early = null;
        this.interactions = null;
        this.rasa = null;
        this.dynamine = null;
        this.efoldmine = null;
        this.z = null;
        this.p = null;

        this.partitionings = new ArrayList<>();
        inSilicoData.forEach((key, value) -> this.partitionings.add(composePartitioning(key, value)));

        this.topology = composeTopology(chain);
    }

    public JsonChain(Chain chain,
                     List<AminoAcid> earlyFoldingResidues,
                     List<Double> dynamineScores,
                     List<Double> efoldmineScores,
                     PartitionedGraph<AminoAcid> experimentalData,
                     Map<String, PartitionedGraph<AminoAcid>> inSilicoData) {
        this.sequence = chain.getAminoAcidSequence();
        this.early = earlyFoldingResidues.stream()
                .map(Group::getResidueIndex)
                .collect(Collectors.toList());
        this.sse = composeSecondaryStructures(chain);
        this.interactions = chain.aminoAcids()
                .map(aminoAcid -> aminoAcid.getFeature(PLIPInteractionContainer.class))
                .map(PLIPInteractionContainer::getInteractions)
                .mapToInt(Collection::size)
                .boxed()
                .collect(Collectors.toList());
        this.energy = chain.aminoAcids()
                .map(aminoAcid -> aminoAcid.getFeature(EnergyProfile.class))
                .mapToDouble(EnergyProfile::getSolvationEnergy)
                .boxed()
                .collect(Collectors.toList());
        this.energy = minMaxNormalize(energy, true);
        this.rasa = chain.aminoAcids()
                .map(aminoAcid -> aminoAcid.getFeature(AccessibleSurfaceArea.class))
                .mapToDouble(AccessibleSurfaceArea::getRelativeAccessibleSurfaceArea)
                .boxed()
                .collect(Collectors.toList());
        this.dynamine = minMaxNormalize(dynamineScores, false);
        this.efoldmine = efoldmineScores;
        PartitionedGraph<AminoAcid> netCartoPartitioning = inSilicoData.get("NetCarto");
        List<double[]> rawPzValues = chain.aminoAcids()
                .map(aminoAcid -> computePz(aminoAcid, netCartoPartitioning))
                .collect(Collectors.toList());
        this.z = rawPzValues.stream()
                .map(values -> values[0])
                .collect(Collectors.toList());
        this.p = rawPzValues.stream()
                .map(values -> values[1])
                .collect(Collectors.toList());

        this.partitionings = new ArrayList<>();
        this.partitionings.add(composePartitioning("HDX-MS Modules", experimentalData));
        inSilicoData.forEach((key, value) -> this.partitionings.add(composePartitioning(key, value)));

        this.topology = null;
    }

    private List<Double> minMaxNormalize(Collection<Double> values, boolean invert) {
        double min = values.stream()
                .mapToDouble(Double::valueOf)
                .min()
                .getAsDouble();
        double max = values.stream()
                .mapToDouble(Double::valueOf)
                .max()
                .getAsDouble();
        return values.stream()
                .mapToDouble(value -> invert ? 1 - (value - min) / (max  - min) :  (value - min) / (max  - min))
                .boxed()
                .collect(Collectors.toList());
    }

    private double[] computePz(AminoAcid aminoAcid, PartitionedGraph<AminoAcid> netCartoPartitioning) {
        Optional<Module<AminoAcid>> s_j_opt = netCartoPartitioning.getModuleOf(aminoAcid);
        if(s_j_opt.isPresent()) {
            Module<AminoAcid> s_j = s_j_opt.get();
            double kappa_i = determineKappaI(aminoAcid, s_j, netCartoPartitioning);
            List<Double> kappa_s_j_values = s_j.getNodes()
                    .stream()
                    .mapToDouble(node -> determineKappaI(node, s_j, netCartoPartitioning))
                    .boxed()
                    .collect(Collectors.toList());
            double kappa_s_j = kappa_s_j_values.stream()
                    .mapToDouble(Double::valueOf)
                    .average()
                    .getAsDouble();
            double sigma_kappa_s_j = Math.sqrt(StatUtils.variance(kappa_s_j_values.stream()
                    .mapToDouble(Double::valueOf)
                    .toArray()));
            double k_i = netCartoPartitioning.getDegreeOf(aminoAcid);

            double z_i = (kappa_i - kappa_s_j) / sigma_kappa_s_j;
            double P_i = 1 - netCartoPartitioning.getModules()
                    .stream()
                    .mapToDouble(module -> {
                        double kappa_is = determineKappaI(aminoAcid, module, netCartoPartitioning);
                        double summand = kappa_is / k_i;
                        return summand * summand;
                    })
                    .sum();

            return new double[] { z_i, P_i };
        } else {
            return new double[] { 0, 0 };
        }
    }

    private double determineKappaI(AminoAcid aminoAcid,
                                   Module<AminoAcid> s_j,
                                   PartitionedGraph<AminoAcid> netCartoPartitioning) {
        return netCartoPartitioning.getEdges()
                .stream()
                .filter(edge -> (edge.getLeft().equals(aminoAcid) && s_j.containsNode(edge.getRight()) ||
                        (edge.getRight().equals(aminoAcid) && s_j.containsNode(edge.getLeft()))))
                .count();
    }

    private JsonPartitioning composePartitioning(String identifier, PartitionedGraph<AminoAcid> data) {
        return new JsonPartitioning(identifier, composeModules(data));
    }

    private List<JsonModule> composeModules(PartitionedGraph<AminoAcid> data) {
        Map<Optional<Module<AminoAcid>>, List<AminoAcid>> map = data.getNodes().stream()
                .collect(Collectors.groupingBy(data::getModuleOf));

        return map.entrySet().stream()
                .filter(entry -> entry.getKey().isPresent())
                .map(this::composeModule)
                .collect(Collectors.toList());
    }

    private JsonModule composeModule(Map.Entry<Optional<Module<AminoAcid>>, List<AminoAcid>> entry) {
        Module<AminoAcid> module = entry.getKey().get();

        List<AminoAcid> aminoAcids = entry.getValue().stream()
                .map(AminoAcid.class::cast)
                .collect(Collectors.toList());

        List<JsonRange> ranges = new ArrayList<>();

        int start = Integer.MIN_VALUE;
        int lastIndex = Integer.MIN_VALUE;

        for(AminoAcid aminoAcid : aminoAcids) {
            int currentIndex = aminoAcid.getResidueIndex() + 1;
            Optional<AminoAcid> previousAminoAcid = aminoAcid.getPreviousAminoAcid();
            if(!previousAminoAcid.isPresent() || !module.containsNode(previousAminoAcid.get())) {
                if (start != Integer.MIN_VALUE) {
                    // terminate last range record if something was observed before
                    ranges.add(new JsonRange(start, lastIndex));
                }
                start = currentIndex;
            }

            lastIndex = currentIndex;
        }
        ranges.add(new JsonRange(start, lastIndex));

        return new JsonModule(module.getIdentifier(), ranges);
    }

    private List<JsonSecondaryStructure> composeSecondaryStructures(Chain chain) {
        List<AminoAcid> aas = chain.aminoAcids().collect(Collectors.toList());
        List<JsonSecondaryStructure> sse = new ArrayList<>();

        String lastSse = "c";
        AminoAcid start = null;
        AminoAcid last = null;
        for(AminoAcid aa : aas) {
            String currentSse = aa.getFeature(GenericSecondaryStructure.class)
                    .getSecondaryStructure()
                    .getReducedRepresentation();
            if(lastSse.equals("c") && !currentSse.equals("c")) {
                // sse observed for the first time
                lastSse = currentSse;
                start = aa;
            } else if(!lastSse.equals(currentSse)) {
                // sse terminated
                if(!start.equals(last)) {
                    sse.add(new JsonSecondaryStructure(start.getResidueIndex() + 1,
                            last.getResidueIndex() + 1,
                            lastSse));
                }
                if(!currentSse.equals("c")) {
                    // sse changed - init new one
                    start = aa;
                }
                lastSse = currentSse;
            }
            last = aa;
        }

        if(!lastSse.equals("c") && !start.equals(last)) {
            sse.add(new JsonSecondaryStructure(start.getResidueIdentifier().getResidueNumber(),
                    last.getResidueIdentifier().getResidueNumber(),
                    lastSse));
        }

        return sse;
    }

    private List<JsonRange> composeTopology(Chain chain) {
        List<AminoAcid> aas = chain.aminoAcids().collect(Collectors.toList());
        List<JsonRange> topology = new ArrayList<>();

        String lastTopology = "o";
        AminoAcid start = null;
        AminoAcid last = null;
        for(AminoAcid aa : aas) {
            String currentTopology = aa.getFeature(Topology.class).isTransmembrane() ? "I" : "o";
            if(lastTopology.equals("o") && !currentTopology.equals("o")) {
                // topology observed for the first time
                lastTopology = currentTopology;
                start = aa;
            } else if(!lastTopology.equals(currentTopology)) {
                // topology terminated
                if(!start.equals(last)) {
                    topology.add(new JsonRange(start.getResidueIndex() + 1,
                            last.getResidueIndex() + 1));
                }
                lastTopology = currentTopology;
            }
            last = aa;
        }

        if(!lastTopology.equals("o") && !start.equals(last)) {
            topology.add(new JsonRange(start.getResidueIdentifier().getResidueNumber(),
                    last.getResidueIdentifier().getResidueNumber()));
        }

        return topology;
    }

    public String getSequence() {
        return sequence;
    }

    public List<Integer> getEarly() {
        return early;
    }

    public List<JsonSecondaryStructure> getSse() {
        return sse;
    }

    public List<Integer> getInteractions() {
        return interactions;
    }

    public List<Double> getEnergy() {
        return energy;
    }

    public List<Double> getRasa() {
        return rasa;
    }

    public List<Double> getDynamine() {
        return dynamine;
    }

    public List<Double> getEfoldmine() {
        return efoldmine;
    }

    public List<Double> getZ() {
        return z;
    }

    public List<Double> getP() {
        return p;
    }

    public List<JsonPartitioning> getPartitionings() {
        return partitionings;
    }

    public List<JsonRange> getTopology() {
        return topology;
    }
}
