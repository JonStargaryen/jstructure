package de.bioforscher.jstructure.membrane.modularity.visualization;

import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.mathematics.graph.PartitionedGraph;
import de.bioforscher.jstructure.mathematics.graph.partitioning.Module;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;

public class JsonChain {
    private final String sequence;
    private final List<JsonSecondaryStructure> sse;
    private final List<JsonPartitioning> partitionings;

    public JsonChain(Chain chain,
                     PartitionedGraph<AminoAcid> experimentalData,
                     Map<String, PartitionedGraph<AminoAcid>> inSilicoData) {
        this.sequence = chain.getAminoAcidSequence();

        this.sse = composeSecondaryStructures(chain);

        this.partitionings = new ArrayList<>();
        this.partitionings.add(composePartitioning("exp", experimentalData));
        inSilicoData.forEach((key, value) -> this.partitionings.add(composePartitioning(key, value)));
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

    public String getSequence() {
        return sequence;
    }

    public List<JsonSecondaryStructure> getSse() {
        return sse;
    }

    public List<JsonPartitioning> getPartitionings() {
        return partitionings;
    }
}
