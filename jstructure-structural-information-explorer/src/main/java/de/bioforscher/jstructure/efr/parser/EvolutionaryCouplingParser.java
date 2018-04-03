package de.bioforscher.jstructure.efr.parser;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.efr.model.HotSpotScoring;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.si.explorer.model.ContactStructuralInformation;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.nodes.Element;
import org.jsoup.select.Elements;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;
import java.util.OptionalDouble;
import java.util.stream.Collectors;

public class EvolutionaryCouplingParser {
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
                                                 Document document) {
        Element table = document.getElementsByTag("table").first();
        contacts.forEach(contact -> parsePlmScore(contact, table));
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

        contact.setPlmScore(plmScore.orElse(0.0));
    }
}
