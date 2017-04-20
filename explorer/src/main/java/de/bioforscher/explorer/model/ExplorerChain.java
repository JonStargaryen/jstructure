package de.bioforscher.explorer.model;

import de.bioforscher.jstructure.mathematics.LinearAlgebraAtom;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.parser.plip.PLIPAnnotator;
import de.bioforscher.jstructure.parser.plip.PLIPInteractionContainer;
import de.bioforscher.jstructure.parser.plip.interaction.PLIPInteraction;
import de.bioforscher.jstructure.parser.uniprot.*;

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
        this.id = chain.getParentProtein().getName().toLowerCase() + "_" + chain.getChainId();
        this.pdb = chain.composePDBRecord();
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
            PLIPInteractionContainer interactions = chain.getParentProtein().getFeature(PLIPInteractionContainer.class, PLIPAnnotator.PLIP_INTERACTIONS);
            LinearAlgebraAtom.Transformation transformation = chain.getFeature(LinearAlgebraAtom.Transformation.class, ModelFactory.TRANSFORMATION);
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

        UniProtAnnotationContainer uniProtAnnotationContainer = chain.getFeature(UniProtAnnotationContainer.class, UniProtAnnotator.UNIPROT_ANNOTATION);
        this.active = uniProtAnnotationContainer.getActiveSites();
        this.disulfide = uniProtAnnotationContainer.getDisulfideBonds();
        this.mutagenesis = uniProtAnnotationContainer.getMutagenesisSites();
        this.ptm = uniProtAnnotationContainer.getAminoAcidModifications();
        this.sse = uniProtAnnotationContainer.getSecondaryStructureElements();
        this.tm = uniProtAnnotationContainer.getTransmembraneRegions();
        this.variants = uniProtAnnotationContainer.getNaturalVariants();
    }

    private List<ExplorerInteraction> convert(List<? extends PLIPInteraction> interactions,
                                              LinearAlgebraAtom.Transformation transformation,
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
