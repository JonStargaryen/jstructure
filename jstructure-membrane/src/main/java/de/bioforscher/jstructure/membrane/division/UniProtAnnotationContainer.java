package de.bioforscher.jstructure.membrane.division;

import org.jsoup.nodes.Document;

import java.util.List;
import java.util.stream.Collectors;

public class UniProtAnnotationContainer {
    private final List<FunctionalSite> activeSites;
//    private final List<TransmembraneRegion> transmembraneRegions;

    UniProtAnnotationContainer(Document document) {
        this.activeSites = document.getElementsByAttribute("type").stream()
                .filter(element -> {
                    String type = element.attr("type");
                    return type.contains("site");
                })
                .map(FunctionalSite::new)
                .collect(Collectors.toList());
//        this.transmembraneRegions = document.getElementsByAttributeValue("type", "transmembrane region").stream()
//                .map(TransmembraneRegion::new)
//                .collect(Collectors.toList());
    }

    public List<FunctionalSite> getActiveSites() {
        return activeSites;
    }

//    public List<TransmembraneRegion> getTransmembraneRegions() {
//        return transmembraneRegions;
//    }
}
