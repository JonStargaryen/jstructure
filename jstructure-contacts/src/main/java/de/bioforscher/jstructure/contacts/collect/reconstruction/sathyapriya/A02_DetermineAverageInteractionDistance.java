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

public class A02_DetermineAverageInteractionDistance {
    private static final Logger logger = LoggerFactory.getLogger(A02_DetermineAverageInteractionDistance.class);

    public static void main(String[] args) {
        Map<? extends Class<? extends PLIPInteraction>, List<PLIPInteraction>> map = ContactsConstants.lines(ContactsConstants.RECONSTRUCTION_DIRECTORY.resolve("ids.list"))
                .flatMap(A02_DetermineAverageInteractionDistance::handlePdbId)
                .collect(Collectors.groupingBy(PLIPInteraction::getClass));

        map.entrySet()
                .stream()
                .map(entry -> entry.getKey().getSimpleName() + " " + StandardFormat.format(entry.getValue()
                        .stream()
                        .mapToDouble(A02_DetermineAverageInteractionDistance::calculateDistance)
                        .average()
                        .getAsDouble()))
                .forEach(System.out::println);

        System.out.println("OverallAverage " + StandardFormat.format(map.values()
                .stream()
                .flatMap(Collection::stream)
                .mapToDouble(A02_DetermineAverageInteractionDistance::calculateDistance)
                .average()
                .getAsDouble()));

        /*
        SaltBridge 6.8865
        WaterBridge 6.5379
        PiCationInteraction 5.8719
        HydrophobicInteraction 5.9347
        PiStacking 4.7513
        HydrogenBond 5.8208
         */

        map.entrySet()
                .stream()
                .map(entry -> entry.getKey().getSimpleName() + " " + StandardFormat.format(entry.getValue()
                        .stream()
                        .mapToDouble(A02_DetermineAverageInteractionDistance::calculateDistance)
                        .max()
                        .getAsDouble()))
                .forEach(System.out::println);

        System.out.println("OverallAverage " + StandardFormat.format(map.values()
                .stream()
                .flatMap(Collection::stream)
                .mapToDouble(A02_DetermineAverageInteractionDistance::calculateDistance)
                .max()
                .getAsDouble()));
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
