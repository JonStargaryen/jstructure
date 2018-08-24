package de.bioforscher.jstructure.efr.parser;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.efr.Start2FoldConstants;
import de.bioforscher.jstructure.efr.model.HotSpotScoring;
import de.bioforscher.jstructure.efr.model.si.ContactStructuralInformation;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.nodes.Element;
import org.jsoup.select.Elements;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.file.Path;
import java.util.Comparator;
import java.util.List;
import java.util.OptionalDouble;
import java.util.OptionalInt;
import java.util.stream.Collectors;

public class EvolutionaryCouplingParser {
    public static void parseHotSpotFile(Chain chain,
                                        Path hotSpotPath) {
        parseHotSpotFile(chain, Start2FoldConstants.newInputStream(hotSpotPath));
    }

    public static void parseHotSpotFile(Chain chain,
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
//                    System.out.println(aminoAcid.getOneLetterCode() + " " + aa);
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

    public static void parsePlmScore(List<ContactStructuralInformation> contacts,
                                     Document document,
                                     int numberOfResidues) {
        Element table = document.getElementsByTag("table").first();
        contacts.forEach(contact -> parsePlmScore(contact, table));

        contacts.sort(Comparator.comparingDouble(ContactStructuralInformation::getCouplingRank));

        double fractionTopScoring02 = 0.2;
        int contactsToSelect02 = (int) (fractionTopScoring02 * numberOfResidues);
        contacts.stream()
                .filter(contactStructuralInformation -> contactStructuralInformation.getCouplingRank() > 0)
                .limit(contactsToSelect02)
                .forEach(ContactStructuralInformation::markAsTopScoringContact02);

        double fractionTopScoring04 = 0.4;
        int contactsToSelect04 = (int) (fractionTopScoring04 * numberOfResidues);
        contacts.stream()
                .filter(contactStructuralInformation -> contactStructuralInformation.getCouplingRank() > 0)
                .limit(contactsToSelect04)
                .forEach(ContactStructuralInformation::markAsTopScoringContact04);

        double fractionTopScoring06 = 0.6;
        int contactsToSelect06 = (int) (fractionTopScoring06 * numberOfResidues);
        contacts.stream()
                .filter(contactStructuralInformation -> contactStructuralInformation.getCouplingRank() > 0)
                .limit(contactsToSelect06)
                .forEach(ContactStructuralInformation::markAsTopScoringContact06);

        double fractionTopScoring08 = 0.8;
        int contactsToSelect08 = (int) (fractionTopScoring08 * numberOfResidues);
        contacts.stream()
                .filter(contactStructuralInformation -> contactStructuralInformation.getCouplingRank() > 0)
                .limit(contactsToSelect08)
                .forEach(ContactStructuralInformation::markAsTopScoringContact08);

        double fractionTopScoring10 = 1.0;
        int contactsToSelect10 = (int) (fractionTopScoring10 * numberOfResidues);
        contacts.stream()
                .filter(contactStructuralInformation -> contactStructuralInformation.getCouplingRank() > 0)
                .limit(contactsToSelect10)
                .forEach(ContactStructuralInformation::markAsTopScoringContact10);

        double fractionTopScoring12 = 1.2;
        int contactsToSelect12 = (int) (fractionTopScoring12 * numberOfResidues);
        contacts.stream()
                .filter(contactStructuralInformation -> contactStructuralInformation.getCouplingRank() > 0)
                .limit(contactsToSelect12)
                .forEach(ContactStructuralInformation::markAsTopScoringContact12);

        double fractionTopScoring14 = 1.4;
        int contactsToSelect14 = (int) (fractionTopScoring14 * numberOfResidues);
        contacts.stream()
                .filter(contactStructuralInformation -> contactStructuralInformation.getCouplingRank() > 0)
                .limit(contactsToSelect14)
                .forEach(ContactStructuralInformation::markAsTopScoringContact14);

        double fractionTopScoring16 = 1.6;
        int contactsToSelect16 = (int) (fractionTopScoring16 * numberOfResidues);
        contacts.stream()
                .filter(contactStructuralInformation -> contactStructuralInformation.getCouplingRank() > 0)
                .limit(contactsToSelect16)
                .forEach(ContactStructuralInformation::markAsTopScoringContact16);
    }

    private static void parsePlmScore(ContactStructuralInformation contact, Element table) {
        String residueIdentifier1 = String.valueOf(contact.getResidueIdentifier1());
        String residueIdentifier2 = String.valueOf(contact.getResidueIdentifier2());

        OptionalDouble plmScore = table.getElementsByTag("tr")
                .stream()
                .skip(1)
                .map(element -> element.getElementsByTag("td"))
                .filter(elements -> (elements.get(2).text().equals(residueIdentifier1) && elements.get(4).text().equals(residueIdentifier2)) ||
                        (elements.get(2).text().equals(residueIdentifier2) && elements.get(4).text().equals(residueIdentifier1)))
                .mapToDouble(elements -> Double.valueOf(elements.get(6).text()))
                .mapToObj(StandardFormat::format)
                .mapToDouble(Double::valueOf)
                .findFirst();
        OptionalInt rank = table.getElementsByTag("tr")
                .stream()
                .skip(1)
                .map(element -> element.getElementsByTag("td"))
                .filter(elements -> (elements.get(2).text().equals(residueIdentifier1) && elements.get(4).text().equals(residueIdentifier2)) ||
                        (elements.get(2).text().equals(residueIdentifier2) && elements.get(4).text().equals(residueIdentifier1)))
                .mapToInt(elements -> Integer.valueOf(elements.get(0).text()))
                .findFirst();

        contact.setPlmScore(plmScore.orElse(0.0));
        contact.setCouplingRank(rank.orElse(-1));
//        contact.setCouplingRank(rank.orElse(Integer.MAX_VALUE));
    }
}
