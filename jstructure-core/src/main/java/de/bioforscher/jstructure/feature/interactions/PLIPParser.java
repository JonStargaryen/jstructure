package de.bioforscher.jstructure.feature.interactions;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import org.jsoup.nodes.Document;
import org.jsoup.nodes.Element;
import org.jsoup.select.Elements;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.lang.reflect.Constructor;
import java.util.ArrayList;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Optional;

/**
 * Parses PLIP files. Especially designed to annotate residue-residue interactions in protein structures.
 * Created by bittrich on 2/9/17.
 */
public class PLIPParser {
    private static final Logger logger = LoggerFactory.getLogger(PLIPParser.class);

    public static List<PLIPInteraction> parse(Chain chain, Document document) {
        List<PLIPInteraction> plipInteractions = new ArrayList<>();
        for(PLIPInteractionType plipInteractionType : PLIPInteractionType.values()) {
            logger.debug("parsing interactions of type {}", plipInteractionType);
            Elements interactionElements = document.getElementsByTag(plipInteractionType.getInteractionTag());
            for(Element interactionElement : interactionElements) {
                Optional<Group> currentGroup = interactionElement.getElementsByTag("resnr")
                        .stream()
                        .map(Element::text)
                        // a safety net for exceedingly big numbers outputted by PLIP
                        .filter(string -> string.length() < 10)
                        .mapToInt(Integer::valueOf)
                        .mapToObj(residueNumber -> chain.select()
                                .residueNumber(residueNumber)
                                .asOptionalGroup())
                        .findFirst()
                        .orElse(Optional.empty());

                if(!currentGroup.isPresent()) {
                    //TODO does a partner in another chain actually make sense?
                    logger.trace("reference to group in different chain or failed to parse line:{}{}", System.lineSeparator(), interactionElement.text());
                    continue;
                }

                try {
                    Constructor<? extends PLIPInteraction> constructor = plipInteractionType.getDescribingClass().getDeclaredConstructor(Group.class, Element.class);
                    constructor.setAccessible(true);
                    PLIPInteraction plipInteraction = constructor.newInstance(currentGroup.get(), interactionElement);
                    plipInteractions.add(plipInteraction);
                } catch (Exception e) {
                    // move to root cause
                    Throwable cause = e;
                    while(cause.getCause() != null) {
                        cause = cause.getCause();
                    }
                    logger.warn("encountered exception during plip parsing: {} - {}", cause.getClass(), cause.getMessage());
                    // interactions whose parsing failed are ignore and do not result in a thrown exception
                    //TODO strict flag
//                    throw new ComputationException(e);
                }
            }
        }

        // merge entries which need merging and remove the merged entries
        List<PLIPInteraction> plipInteractionsToRemove = new ArrayList<>();
        for(PLIPInteraction plipInteraction : plipInteractions) {
            Group partner1 = plipInteraction.getPartner1();
            Group partner2 = plipInteraction.getPartner2();

            if(plipInteraction instanceof PiCationInteraction || plipInteraction instanceof PiStacking || plipInteraction instanceof SaltBridge) {
                // keep only those interactions where the first resNum is smaller
                if(partner1.getResidueNumber().getResidueNumber() > partner2.getResidueNumber().getResidueNumber()) {
                    plipInteractionsToRemove.add(plipInteraction);
                    continue;
                }

                try {
                    if (plipInteraction instanceof PiCationInteraction) {
                        PiCationInteraction otherHalfOfInteraction = plipInteractions.stream()
                                .filter(PiCationInteraction.class::isInstance)
                                .map(PiCationInteraction.class::cast)
                                .filter(piCationInteraction -> piCationInteraction.getPartner1().equals(partner2) && piCationInteraction.getPartner2().equals(partner1))
                                .findFirst()
                                .orElseThrow(NoSuchElementException::new);
                        ((PiCationInteraction) plipInteraction).getAtoms1().addAll(otherHalfOfInteraction.getAtoms2());
                    } else if (plipInteraction instanceof PiStacking) {
                        PiStacking otherHalfOfInteraction = plipInteractions.stream()
                                .filter(PiStacking.class::isInstance)
                                .map(PiStacking.class::cast)
                                .filter(piStacking -> piStacking.getPartner1().equals(partner2) && piStacking.getPartner2().equals(partner1))
                                .findFirst()
                                .orElseThrow(NoSuchElementException::new);
                        ((PiStacking) plipInteraction).getAtoms1().addAll(otherHalfOfInteraction.getAtoms2());
                    } else {
                        SaltBridge otherHalfOfInteraction = plipInteractions.stream()
                                .filter(SaltBridge.class::isInstance)
                                .map(SaltBridge.class::cast)
                                .filter(saltBridge -> saltBridge.getPartner1().equals(partner2) && saltBridge.getPartner2().equals(partner1))
                                .findFirst()
                                .orElseThrow(NoSuchElementException::new);
                        ((SaltBridge) plipInteraction).getAtoms1().addAll(otherHalfOfInteraction.getAtoms2());
                    }
                } catch (NoSuchElementException e) {
                    logger.debug("could not find other half of {} {} {}", plipInteraction.getClass().getSimpleName(), plipInteraction.getPartner1().getIdentifier(), plipInteraction.getPartner2().getIdentifier());
                }
            }
        }

        plipInteractions.removeAll(plipInteractionsToRemove);

        return plipInteractions;
    }
}
