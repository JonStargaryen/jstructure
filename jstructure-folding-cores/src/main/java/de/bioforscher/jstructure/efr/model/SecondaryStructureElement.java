package de.bioforscher.jstructure.efr.model;

import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

public class SecondaryStructureElement {
    private final int start;
    private final int end;
    private final String type;

    public SecondaryStructureElement(int start, int end, String type) {
        this.start = start;
        this.end = end;
        this.type = type;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    public String getType() {
        return type;
    }

    public static List<SecondaryStructureElement> of(Chain chain) {
        List<SecondaryStructureElement> secondaryStructureElements = new ArrayList<>();
        int start = -1;
        int end = -1;
        String currentType = "c";
        List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());

        for(AminoAcid aminoAcid : aminoAcids) {
            String type = aminoAcid.getFeature(GenericSecondaryStructure.class)
                    .getSecondaryStructure()
                    .getReducedRepresentation();
            if(!currentType.equals("c") && !currentType.equals(type)) {
                secondaryStructureElements.add(new SecondaryStructureElement(start,
                        end,
                        currentType));
                start = -1;
                end = -1;
            }

            currentType = type;
            if(!currentType.equals("c")) {
                int resnum = aminoAcid.getResidueIndex();
//                int resnum = aminoAcid.getResidueIdentifier().getResidueNumber();
                if(start == -1) {
                    start = resnum;
                }
                end = resnum;
            }
        }

        if(!currentType.equals("c") && start != -1) {
            secondaryStructureElements.add(new SecondaryStructureElement(start,
                    end,
                    currentType));
        }
        return secondaryStructureElements;
    }
}
