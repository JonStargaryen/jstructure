package de.bioforscher.jstructure.contacts.collect.reconstruction.sathyapriya;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.contacts.ContactsConstants;
import de.bioforscher.jstructure.feature.interactions.PLIPInteraction;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class A02_DetermineInteractionDistances {
    private static final Logger logger = LoggerFactory.getLogger(A02_DetermineInteractionDistances.class);

    public static void main(String[] args) {
        Map<? extends Class<? extends PLIPInteraction>, List<PLIPInteraction>> map = ContactsConstants.lines(ContactsConstants.RECONSTRUCTION_DIRECTORY.resolve("ids.list"))
                .flatMap(A02_DetermineInteractionDistances::handlePdbId)
                .collect(Collectors.groupingBy(PLIPInteraction::getClass));

        System.out.println("average");
        map.entrySet()
                .stream()
                .map(entry -> entry.getKey().getSimpleName() + " " + StandardFormat.format(entry.getValue()
                        .stream()
                        .mapToDouble(A02_DetermineInteractionDistances::calculateDistance)
                        .average()
                        .getAsDouble()))
                .forEach(System.out::println);

        System.out.println("Generic " + StandardFormat.format(map.values()
                .stream()
                .flatMap(Collection::stream)
                .mapToDouble(A02_DetermineInteractionDistances::calculateDistance)
                .average()
                .getAsDouble()));

        /*
        Averages:
        SaltBridge 6.8865
        WaterBridge 6.5379
        HydrophobicInteraction 5.9347
        PiStacking 4.7513
        PiCationInteraction 5.8719
        HydrogenBond 5.8208
        Generic 5.9308

        performance:
        - degenerates performance compared to max values
         */

        System.out.println();
        System.out.println("max");
        map.entrySet()
                .stream()
                .map(entry -> entry.getKey().getSimpleName() + " " + StandardFormat.format(entry.getValue()
                        .stream()
                        .mapToDouble(A02_DetermineInteractionDistances::calculateDistance)
                        .max()
                        .getAsDouble()))
                .forEach(System.out::println);

        System.out.println("Generic " + StandardFormat.format(map.values()
                .stream()
                .flatMap(Collection::stream)
                .mapToDouble(A02_DetermineInteractionDistances::calculateDistance)
                .max()
                .getAsDouble()));

        /*
        Maximum values:
        SaltBridge 10.8307
        WaterBridge 12.1106
        HydrophobicInteraction 11.9862
        PiStacking 5.0976
        PiCationInteraction 7.393
        HydrogenBond 11.5943
        Generic 12.1106

        performance:
        - less contacts, lower performance
         */

        System.out.println();
        System.out.println("backbone interactions");
        map.entrySet()
                .stream()
                .map(entry -> entry.getKey().getSimpleName() + " " + StandardFormat.format(entry.getValue()
                        .stream()
                        .filter(PLIPInteraction::isBackboneInteraction)
                        .mapToDouble(A02_DetermineInteractionDistances::calculateDistance)
                        .max()
                        .orElse(-1)))
                .forEach(System.out::println);

        System.out.println("Generic " + StandardFormat.format(map.values()
                .stream()
                .flatMap(Collection::stream)
                .mapToDouble(A02_DetermineInteractionDistances::calculateDistance)
                .max()
                .getAsDouble()));

        /*
        Backbone interactions:
        SaltBridge -1
        WaterBridge 8.8893
        HydrophobicInteraction -1
        PiStacking -1
        PiCationInteraction -1
        HydrogenBond 10.0201
        Generic 12.1106
        - not applicable as additional constraints
        */
    }

    private static double calculateDistance(PLIPInteraction plipInteraction) {
        AminoAcid partner1 = (AminoAcid) plipInteraction.getPartner1();
        AminoAcid partner2 = (AminoAcid) plipInteraction.getPartner2();

        return ContactsConstants.getBetaCarbon(partner1).calculate().distance(ContactsConstants.getBetaCarbon(partner2));
    }

    private static Stream<PLIPInteraction> handlePdbId(String pdbId) {
        Structure structure = StructureParser.source(pdbId).parse();
        Chain chain = structure.chains().findFirst().get();
        logger.info("processing {} with {} chains",
                pdbId,
                structure.chains().count());

        PLIPInteractionContainer interactionContainer = chain.getFeature(PLIPInteractionContainer.class);
        return interactionContainer.getInteractions().stream();
    }
}
