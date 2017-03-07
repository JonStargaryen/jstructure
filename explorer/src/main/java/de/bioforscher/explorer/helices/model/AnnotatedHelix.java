package de.bioforscher.explorer.helices.model;

/**
 * Represent one helix from the AHAH data set.
 * Created by bittrich on 3/7/17.
 */
public class AnnotatedHelix {
    private String helixId, pdbId, chainId;
    private double kinked, curved, straight;
    private int start, end, total;
    private HelixClassification classification;
    private HelixProtein protein;

    public AnnotatedHelix(String helixId, String pdbId, String chainId, int start, int end, double kinked, double curved, double straight, int total, HelixClassification classification, HelixProtein protein) {
        this.helixId = helixId;
        this.pdbId = pdbId;
        this.chainId = chainId;
        this.start = start;
        this.end = end;
        this.kinked = kinked;
        this.curved = curved;
        this.straight = straight;
        this.total = total;
        this.classification = classification;
        this.protein = protein;
    }

    public AnnotatedHelix() {

    }

    public String getHelixId() {
        return helixId;
    }

    public String getPdbId() {
        return pdbId;
    }

    public String getChainId() {
        return chainId;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    public double getKinked() {
        return kinked;
    }

    public double getCurved() {
        return curved;
    }

    public double getStraight() {
        return straight;
    }

    public int getTotal() {
        return total;
    }

    public HelixClassification getClassification() {
        return classification;
    }

    public HelixProtein getProtein() {
        return protein;
    }
}
