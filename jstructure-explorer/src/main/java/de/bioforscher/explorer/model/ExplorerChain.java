package de.bioforscher.explorer.model;

import de.bioforscher.jstructure.feature.interactions.PLIPInteraction;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.uniprot.*;
import de.bioforscher.jstructure.mathematics.Transformation;
import de.bioforscher.jstructure.model.structure.Chain;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * One protein chain in the explorer.
 * Created by bittrich on 4/6/17.
 */
public class ExplorerChain {
    private String id, pdb, rep, sequence;
    private List<ExplorerLigand> ligands;
    private boolean isRep, isPlip;
    /* 'traditional' data */
    private List<ExplorerGroup> groups;
    /* interactions */
    private List<ExplorerInteraction> halogenBonds;
    private List<ExplorerInteraction> hydrogenBonds;
    private List<ExplorerInteraction> metalComplexes;
    private List<ExplorerInteraction> piCationInteractions;
    private List<ExplorerInteraction> piStackings;
    private List<ExplorerInteraction> saltBridges;
    private List<ExplorerInteraction> waterBridges;
    /* UniProt data */
    private List<UniProtActiveSite> active;
    private List<UniProtDisulfideBond> disulfide;
    private List<UniProtMutagenesisSite> mutagenesis;
    private List<UniProtAminoAcidModification> ptm;
    private List<UniProtSecondaryStructureElement> sse;
    private List<UniProtTransmembraneRegion> tm;
    private List<UniProtNaturalVariant> variants;

    public ExplorerChain() {
    }

    public ExplorerChain(Chain chain, String representativeId) {
        this.id = chain.getChainId().getFullName();
        this.pdb = chain.getPdbRepresentation();
        this.rep = representativeId;
        this.sequence = chain.getAminoAcidSequence();
        this.isRep = representativeId.equals(this.id);

        this.ligands = chain.select()
                .hetatms()
                .asFilteredGroups()
                .map(ExplorerLigand::new)
                .collect(Collectors.toList());

        this.groups = chain.aminoAcids()
                .map(ExplorerGroup::new)
                .collect(Collectors.toList());

        try {
            PLIPInteractionContainer interactions = chain.getFeatureContainer().getFeature(PLIPInteractionContainer.class);
            Transformation transformation = chain.getFeatureContainer().getFeature(Transformation.class);
            this.halogenBonds = convert(interactions.getHalogenBonds(), transformation, chain);
            this.hydrogenBonds = convert(interactions.getHydrogenBonds(), transformation, chain);
            this.metalComplexes = convert(interactions.getMetalComplexes(), transformation, chain);
            this.piCationInteractions = convert(interactions.getPiCationInteractions(), transformation, chain);
            this.piStackings = convert(interactions.getPiStackings(), transformation, chain);
            this.saltBridges = convert(interactions.getSaltBridges(), transformation, chain);
            this.waterBridges = convert(interactions.getWaterBridges(), transformation, chain);
            this.isPlip = true;
        } catch (NullPointerException e) {
            this.halogenBonds = new ArrayList<>();
            this.hydrogenBonds = new ArrayList<>();
            this.metalComplexes = new ArrayList<>();
            this.piCationInteractions = new ArrayList<>();
            this.piStackings = new ArrayList<>();
            this.saltBridges = new ArrayList<>();
            this.waterBridges = new ArrayList<>();
            this.isPlip = false;
        }

        try {
            UniProtAnnotationContainer uniProtAnnotationContainer = chain.getFeatureContainer().getFeature(UniProtAnnotationContainer.class);
            this.active = uniProtAnnotationContainer.getActiveSites();
            this.disulfide = uniProtAnnotationContainer.getDisulfideBonds();
            this.mutagenesis = uniProtAnnotationContainer.getMutagenesisSites();
            this.ptm = uniProtAnnotationContainer.getAminoAcidModifications();
            this.sse = uniProtAnnotationContainer.getSecondaryStructureElements();
            this.tm = uniProtAnnotationContainer.getTransmembraneRegions();
            this.variants = uniProtAnnotationContainer.getNaturalVariants();
        } catch (NullPointerException e) {
            this.active = new ArrayList<>();
            this.disulfide = new ArrayList<>();
            this.mutagenesis = new ArrayList<>();
            this.ptm = new ArrayList<>();
            this.sse = new ArrayList<>();
            this.tm = new ArrayList<>();
            this.variants = new ArrayList<>();
        }
    }

    private List<ExplorerInteraction> convert(List<? extends PLIPInteraction> interactions,
                                              Transformation transformation,
                                              Chain chain) {
        return interactions.stream()
                .filter(interaction -> interaction.getPartner1().getParentChain().getChainId().equals(chain.getChainId()) && interaction.getPartner2().getParentChain().getChainId().equals(chain.getChainId()))
                .map(interaction -> new ExplorerInteraction(interaction, transformation))
                .collect(Collectors.toList());
    }

    public String getId() {
        return id;
    }

    public String getPdb() {
        return pdb;
    }

    public String getRep() {
        return rep;
    }

    public String getSequence() {
        return sequence;
    }

    public List<ExplorerLigand> getLigands() {
        return ligands;
    }

    public boolean isRep() {
        return isRep;
    }

    public boolean isPlip() {
        return isPlip;
    }

    public List<ExplorerInteraction> getHalogenBonds() {
        return halogenBonds;
    }

    public List<ExplorerInteraction> getHydrogenBonds() {
        return hydrogenBonds;
    }

    public List<ExplorerInteraction> getMetalComplexes() {
        return metalComplexes;
    }

    public List<ExplorerInteraction> getPiCationInteractions() {
        return piCationInteractions;
    }

    public List<ExplorerInteraction> getPiStackings() {
        return piStackings;
    }

    public List<ExplorerInteraction> getSaltBridges() {
        return saltBridges;
    }

    public List<ExplorerInteraction> getWaterBridges() {
        return waterBridges;
    }

    public List<ExplorerGroup> getGroups() {
        return groups;
    }

    public List<UniProtActiveSite> getActive() {
        return active;
    }

    public List<UniProtDisulfideBond> getDisulfide() {
        return disulfide;
    }

    public List<UniProtMutagenesisSite> getMutagenesis() {
        return mutagenesis;
    }

    public List<UniProtAminoAcidModification> getPtm() {
        return ptm;
    }

    public List<UniProtSecondaryStructureElement> getSse() {
        return sse;
    }

    public List<UniProtTransmembraneRegion> getTm() {
        return tm;
    }

    public List<UniProtNaturalVariant> getVariants() {
        return variants;
    }
}
