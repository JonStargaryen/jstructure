package de.bioforscher.jstructure.feature.interaction;

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
import java.util.Optional;

public class PLIPLigandParser {
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
                    logger.info("reference to group in different chain or failed to parse line:{}{}", System.lineSeparator(), interactionElement.text());
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

        return plipInteractions;
    }
}