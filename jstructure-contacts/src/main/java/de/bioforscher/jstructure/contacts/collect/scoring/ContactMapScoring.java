package de.bioforscher.jstructure.contacts.collect.scoring;

import de.bioforscher.jstructure.mathematics.SetOperations;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.util.List;

public class ContactMapScoring {
    private final ConfusionMatrix conventionalConfusionMatrix;
    private final ConfusionMatrix plipConfusionMatrix;

    public ContactMapScoring(List<AminoAcid> aminoAcids,
                             ContactMap conventionalReferenceContactMap,
                             ContactMap conventionalModelContactMap,
                             ContactMap plipReferenceContactMap,
                             ContactMap plipModelContactMap) {
        this.conventionalConfusionMatrix = new ConfusionMatrix(aminoAcids, conventionalReferenceContactMap, conventionalModelContactMap);
        this.plipConfusionMatrix = new ConfusionMatrix(aminoAcids, plipReferenceContactMap, plipModelContactMap);
    }

    public ConfusionMatrix getConventionalConfusionMatrix() {
        return conventionalConfusionMatrix;
    }

    public ConfusionMatrix getPlipConfusionMatrix() {
        return plipConfusionMatrix;
    }

    public static class ConfusionMatrix {
        private final double truePositive;
        private final double trueNegative;
        private final double falsePositive;
        private final double falseNegative;

        ConfusionMatrix(List<AminoAcid> aminoAcids,
                        ContactMap reference,
                        ContactMap model) {
            double tp;
            double tn;
            double fp;
            double fn;
            try {
                tp = (int) reference.getContacts()
                        .stream()
                        .filter(model::hasContact)
                        .count();
                tn = (int) SetOperations.uniquePairsOf(aminoAcids)
                        .filter(pair -> !reference.hasContact(pair))
                        .filter(pair -> !model.hasContact(pair))
                        .count();
                fp = (int) SetOperations.uniquePairsOf(aminoAcids)
                        .filter(pair -> !reference.hasContact(pair))
                        .filter(model::hasContact)
                        .count();
                fn = (int) SetOperations.uniquePairsOf(aminoAcids)
                        .filter(reference::hasContact)
                        .filter(pair -> !model.hasContact(pair))
                        .count();
            } catch (Exception e) {
                tp = Double.NaN;
                tn = Double.NaN;
                fp = Double.NaN;
                fn= Double.NaN;
            }
            this.truePositive = tp;
            this.trueNegative = tn;
            this.falsePositive = fp;
            this.falseNegative = fn;
        }

        public double getTruePositive() {
            return truePositive;
        }

        public double getTrueNegative() {
            return trueNegative;
        }

        public double getFalsePositive() {
            return falsePositive;
        }

        public double getFalseNegative() {
            return falseNegative;
        }

        public double getSensitivity() {
            return truePositive / (truePositive + falseNegative);
        }

        public double getSpecificity() {
            return trueNegative / (trueNegative + falsePositive);
        }

        public double getAccuracy() {
            return (truePositive + trueNegative) / (truePositive + trueNegative + falsePositive + falseNegative);
        }

        public double getF1Score() {
            double tp2 = 2 * truePositive;
            return tp2 / (tp2 + falsePositive + falseNegative);
        }

        public double getMCC() {
            return (truePositive * trueNegative + falsePositive * falseNegative) /
                    Math.sqrt((truePositive + falsePositive) * (truePositive + falseNegative) * (trueNegative + falsePositive) * (trueNegative + falseNegative));
        }
    }
}
