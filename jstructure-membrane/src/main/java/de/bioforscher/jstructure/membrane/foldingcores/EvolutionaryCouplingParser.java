package de.bioforscher.jstructure.membrane.foldingcores;

import de.bioforscher.jstructure.contacts.ContactsConstants;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.select.Elements;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.file.Path;
import java.util.List;
import java.util.stream.Collectors;

public class EvolutionaryCouplingParser {
    static void parseHotSpotFile(Chain chain,
                                 InputStream inputStream) {
        Document hotSpotDocument = Jsoup.parse(new BufferedReader(new InputStreamReader(inputStream))
                .lines()
                .collect(Collectors.joining(System.lineSeparator())));
        List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());

        hotSpotDocument.getElementsByTag("tr").stream()
                .skip(1)
                .forEach(element -> {
                    Elements tds = element.getElementsByTag("td");
                    int residueNumber = Integer.valueOf(tds.get(0).text());
                    String aa = tds.get(1).text();
                    int ecCount = Integer.valueOf(tds.get(2).text());
                    double cumStrength = Double.valueOf(tds.get(3).text());
                    double ecStrength = Double.valueOf(tds.get(4).text());
                    int conservation = Integer.valueOf(tds.get(5).text());

                    AminoAcid aminoAcid = aminoAcids.get(residueNumber - 1);
                    aminoAcid.getFeatureContainer().addFeature(new HotSpotScoring(ecCount,
                            cumStrength,
                            ecStrength,
                            conservation));
                });

        // assign baseline
        aminoAcids.stream()
                .map(AminoAcid::getFeatureContainer)
                .filter(featureContainer -> !featureContainer.getFeatureOptional(HotSpotScoring.class).isPresent())
                .forEach(featureContainer -> featureContainer.addFeature(new HotSpotScoring()));
    }

    public static void parseHotSpotFile(Chain chain,
                                        Path hotSpotPath) {
        parseHotSpotFile(chain, ContactsConstants.newInputStream(hotSpotPath));
    }

    public static class HotSpotScoring extends FeatureContainerEntry {
        private final int ecCount;
        private final double cumStrength;
        private final double ecStrength;
        private final int conservation;

        HotSpotScoring(int ecCount, double cumStrength, double ecStrength, int conservation) {
            super(null);
            this.ecCount = ecCount;
            this.cumStrength = cumStrength;
            this.ecStrength = ecStrength;
            this.conservation = conservation;
        }

        HotSpotScoring() {
            this(0, 0, 0, 0);
        }

        public int getEcCount() {
            return ecCount;
        }

        public double getCumStrength() {
            return cumStrength;
        }

        public double getEcStrength() {
            return ecStrength;
        }

        public int getConservation() {
            return conservation;
        }
    }
}
