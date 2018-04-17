package de.bioforscher.jstructure.feature.sse;

import de.bioforscher.jstructure.feature.sse.dssp.DictionaryOfProteinSecondaryStructure;
import de.bioforscher.jstructure.model.feature.ComputationException;
import de.bioforscher.jstructure.model.feature.DefaultFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.util.stream.Collectors;

/**
 * The wrapper for secondary structure elements.
 * Created by bittrich on 6/28/17.
 */
@DefaultFeatureProvider(DictionaryOfProteinSecondaryStructure.class)
public class GenericSecondaryStructure extends FeatureContainerEntry {
    protected SecondaryStructureType secondaryStructure;

    public GenericSecondaryStructure(FeatureProvider featureProvider, SecondaryStructureType secondaryStructure) {
        super(featureProvider);
        this.secondaryStructure = secondaryStructure;
    }

    public SecondaryStructureType getSecondaryStructure() {
        return secondaryStructure;
    }

    public void setSecondaryStructure(SecondaryStructureType secondaryStructure) {
        this.secondaryStructure = secondaryStructure;
    }

    public SecondaryStructureElement getSurroundingSecondaryStructureElement(AminoAcid aminoAcid) {
        return new SecondaryStructureElement(aminoAcid);
    }

    public static class SecondaryStructureElement {
        private final String reducedType;
        private final int size;
        private final int start;
        private final int end;
        private final int terminusDistance;

        SecondaryStructureElement(AminoAcid aminoAcid) {
            Chain chain = aminoAcid.getParentChain();
            String sseString = chain.aminoAcids()
                    .map(aa -> aa.getFeature(GenericSecondaryStructure.class).getSecondaryStructure().getReducedRepresentation())
                    .collect(Collectors.joining());

            // safety net to ensure that numbering matches
            if(sseString.length() != chain.aminoAcids().count()) {
                throw new ComputationException("group numbering was compromised - do not skip ligand parsing when using Group#getResidueIndex");
            }

            reducedType = aminoAcid.getFeature(GenericSecondaryStructure.class).getSecondaryStructure().getReducedRepresentation();

            char aaSseChar = reducedType.charAt(0);
            int index = aminoAcid.getResidueIndex();
            int left = 0;
            int right = 0;
            try {
                for (int i = index - 1; index >= 0; i--) {
                    if (i < 0 || sseString.charAt(i) != aaSseChar) {
                        break;
                    }
                    left++;
                }
                for (int i = index + 1; index < sseString.length(); i++) {
                    if (i >= sseString.length() || sseString.charAt(i) != aaSseChar) {
                        break;
                    }
                    right++;
                }
                this.start = index - left;
                this.end = index + right;
                this.size = left + right + 1;
                this.terminusDistance = Math.min(index - start, end - index);
            } catch (StringIndexOutOfBoundsException e) {
                throw new ComputationException("group numbering was compromised - do not skip ligand parsing when using Group#getResidueIndex", e);
            }
        }

        public String getReducedType() {
            return reducedType;
        }

        public int getSize() {
            return size;
        }

        public int getTerminusDistance() {
            return terminusDistance;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            SecondaryStructureElement that = (SecondaryStructureElement) o;

            if (size != that.size) return false;
            if (start != that.start) return false;
            if (end != that.end) return false;
            return reducedType != null ? reducedType.equals(that.reducedType) : that.reducedType == null;
        }

        @Override
        public int hashCode() {
            int result = reducedType != null ? reducedType.hashCode() : 0;
            result = 31 * result + size;
            result = 31 * result + start;
            result = 31 * result + end;
            return result;
        }

        @Override
        public String toString() {
            return "SecondaryStructureElement{" +
                    "reducedType='" + reducedType + '\'' +
                    ", size=" + size +
                    ", start=" + start +
                    ", end=" + end +
                    '}';
        }
    }
}
