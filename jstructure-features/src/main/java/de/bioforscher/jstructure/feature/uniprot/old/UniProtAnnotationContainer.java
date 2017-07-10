package de.bioforscher.jstructure.feature.uniprot.old;

import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.structure.Chain;
import org.jsoup.nodes.Document;

import java.util.List;
import java.util.stream.Collectors;

/**
 * Gathers all UniProt information of a {@link Chain}.
 * Created by bittrich on 3/2/17.
 */
@Deprecated
public class UniProtAnnotationContainer extends FeatureContainerEntry {
    private String uniProtId;
    private String uniProtSequence;
    private List<UniProtAminoAcidModification> aminoAcidModifications;
    private List<UniProtActiveSite> activeSites;
    private List<UniProtDisulfideBond> disulfideBonds;
    private List<UniProtMutagenesisSite> mutagenesisSites;
    private List<UniProtNaturalVariant> naturalVariants;
    private List<UniProtReference> references;
    private List<UniProtSecondaryStructureElement> secondaryStructureElements;
    private List<UniProtTransmembraneRegion> transmembraneRegions;
    //TODO motifs?

    UniProtAnnotationContainer(AbstractFeatureProvider featureProvider, String uniProtId, Document describingDocument) {
        super(featureProvider);
        this.uniProtId = uniProtId;
        this.uniProtSequence = describingDocument.getElementsByTag("sequence").text().replaceAll("\\s+", "");
        this.aminoAcidModifications = describingDocument.getElementsByTag("feature").stream()
                //TODO register further modifications
                .filter(element -> element.hasAttr("type"))
                .filter(element -> {
                    String type = element.attr("type");
                    return type.equals("lipid moiety-binding region") ||
                            type.equals("glycosylation site") ||
                            type.equals("modified residue");
                })
                .map(UniProtAminoAcidModification::new)
                .collect(Collectors.toList());
        this.activeSites = describingDocument.getElementsByAttribute("type").stream()
                .filter(element -> {
                    String type = element.attr("type");
                    //TODO register further active sites
                    return type.equals("active site") || type.contains("binding site");
                })
                .map(UniProtActiveSite::new)
                .collect(Collectors.toList());
        this.disulfideBonds = describingDocument.getElementsByAttributeValue("type", "disulfide bond").stream()
                // there are also interchain bonds annotated - ignoring them for now
                .filter(element -> !element.hasAttr("description") && element.getElementsByTag("location").size() > 1)
                .map(UniProtDisulfideBond::new)
                .collect(Collectors.toList());
        this.mutagenesisSites = describingDocument.getElementsByAttributeValue("type", "mutagenesis site").stream()
                .map(UniProtMutagenesisSite::new)
                .collect(Collectors.toList());
        this.naturalVariants = describingDocument.getElementsByAttributeValue("type", "sequence variant").stream()
                .map(UniProtNaturalVariant::new)
                .collect(Collectors.toList());
        this.references = describingDocument.getElementsByTag("reference").stream()
                // some reference are no real papers, keep them to keep evidence mapping intact
//                .filter(element -> !element.getElementsByTag("citation").first().attr("type").equals("submission"))
                .map(UniProtReference::new)
                .collect(Collectors.toList());
        this.secondaryStructureElements = describingDocument.getElementsByTag("feature").stream()
                .filter(element -> element.hasAttr("type"))
                .filter(element -> element.attr("type").equals("strand") || element.attr("type").equals("helix"))
                .map(UniProtSecondaryStructureElement::new)
                .collect(Collectors.toList());
        this.transmembraneRegions = describingDocument.getElementsByAttributeValue("type", "transmembrane region").stream()
                .map(UniProtTransmembraneRegion::new)
                .collect(Collectors.toList());
    }

    public String getUniProtId() {
        return uniProtId;
    }

    public String getUniProtSequence() {
        return uniProtSequence;
    }

    public List<UniProtAminoAcidModification> getAminoAcidModifications() {
        return aminoAcidModifications;
    }

    public List<UniProtActiveSite> getActiveSites() {
        return activeSites;
    }

    public List<UniProtDisulfideBond> getDisulfideBonds() {
        return disulfideBonds;
    }

    public List<UniProtMutagenesisSite> getMutagenesisSites() {
        return mutagenesisSites;
    }

    public List<UniProtNaturalVariant> getNaturalVariants() {
        return naturalVariants;
    }

    public List<UniProtReference> getReferences() {
        return references;
    }

    public List<UniProtSecondaryStructureElement> getSecondaryStructureElements() {
        return secondaryStructureElements;
    }

    public List<UniProtTransmembraneRegion> getTransmembraneRegions() {
        return transmembraneRegions;
    }
}
