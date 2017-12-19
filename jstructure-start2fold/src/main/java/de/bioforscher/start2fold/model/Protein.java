package de.bioforscher.start2fold.model;

import org.jsoup.nodes.Element;
import org.jsoup.select.Elements;

public class Protein {
    private final String name;
    private final String entryId;
    private final String uniProtId;
    private final String pdbId;
    private final int numberOfResidues;
    private final String fragment;
    private final String secondaryStructure;

    public Protein(Element element) {
        this(element.getElementsByTag("td"));
    }


    public Protein(Elements elements) {
        this(elements.get(0).text(),
                elements.get(1).text(),
                elements.get(2).text(),
                elements.get(3).text(),
                Integer.valueOf(elements.get(4).text()),
                elements.get(5).text(),
                elements.get(6).text());
    }

    public Protein(String name,
                   String entryId,
                   String uniProtId,
                   String pdbId,
                   int numberOfResidues,
                   String fragment,
                   String secondaryStructure) {
        this.name = name;
        this.entryId = entryId;
        this.uniProtId = uniProtId;
        this.pdbId = pdbId;
        this.numberOfResidues = numberOfResidues;
        this.fragment = fragment;
        this.secondaryStructure = secondaryStructure;
    }

    public String getName() {
        return name;
    }

    public String getEntryId() {
        return entryId;
    }

    public String getUniProtId() {
        return uniProtId;
    }

    public String getPdbId() {
        return pdbId;
    }

    public int getNumberOfResidues() {
        return numberOfResidues;
    }

    public String getFragment() {
        return fragment;
    }

    public String getSecondaryStructure() {
        return secondaryStructure;
    }

    @Override
    public String toString() {
        return "Protein{" +
                "name='" + name + '\'' +
                ", entryId='" + entryId + '\'' +
                ", uniProtId='" + uniProtId + '\'' +
                ", pdbId='" + pdbId + '\'' +
                ", numberOfResidues=" + numberOfResidues +
                ", fragment='" + fragment + '\'' +
                ", secondaryStructure='" + secondaryStructure + '\'' +
                '}';
    }
}
