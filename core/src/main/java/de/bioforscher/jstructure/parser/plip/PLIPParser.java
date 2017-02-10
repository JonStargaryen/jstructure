package de.bioforscher.jstructure.parser.plip;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.nodes.Element;
import org.jsoup.select.Elements;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

/**
 * Parses PLIP files. Especially designed to annotate residue-residue interactions in protein structures.
 * Created by bittrich on 2/9/17.
 */
class PLIPParser {
    private static final Logger logger = LoggerFactory.getLogger(PLIPParser.class);
    private static final String RESNR_TAG = "resnr";

    static List<PLIPInteraction> parse(Chain chain, String documentDescribingChain) {
        Document document = Jsoup.parse(documentDescribingChain);

        List<PLIPInteraction> plipInteractions = new ArrayList<>();
        for(PLIPInteractionType plipInteractionType : PLIPInteractionType.values()) {
            logger.debug("parsing interactions of type {}", plipInteractionType);
            Elements interactionElements = document.getElementsByTag(plipInteractionType.getInteractionTag());
            for(Element interactionElement : interactionElements) {
                Optional<Group> currentGroup = interactionElement.getElementsByTag(RESNR_TAG)
                        .stream()
                        .map(Element::text)
                        // a safety net for exceedingly big numbers outputted by PLIP
                        .filter(string -> string.length() < 10)
                        .map(Integer::valueOf)
                        .map(residueNumber -> Selection.on(chain)
                                .residueNumber(residueNumber)
                                .asOptionalGroup())
                        .findFirst()
                        .orElse(Optional.empty());

                if(!currentGroup.isPresent()) {
                    logger.trace("reference to group in different chain or failed to parse line:" + System.lineSeparator() + interactionElement.text());
                    continue;
                }

                // list of all atoms annotated as interacting with the current group by the current interaction type
                List<Atom> interactingAtoms = new ArrayList<>();

                for(String atomTag : plipInteractionType.getAtomTags()) {
                    Elements interactingAtomTags = interactionElement.getElementsByTag(atomTag);
                    // based on type decide what to do
                    switch (plipInteractionType) {
                        case SALT_BRIDGE: case PI_STACKING: case PI_CATION:
                            Element saltBridgeAtoms = interactingAtomTags.first();
                            interactingAtomTags = saltBridgeAtoms.getElementsByTag("idx");
                        default:
                            interactingAtomTags.stream()
                                    .map(Element::text)
                                    .map(Integer::valueOf)
                                    .map(atomNumber -> Selection.on(chain)
                                            .pdbSerial(atomNumber)
                                            .asOptionalAtom())
                                    .filter(Optional::isPresent)
                                    .map(Optional::get)
                                    .forEach(interactingAtoms::add);
                    }
                }

                // map atoms to their parent group
                List<Group> interactingGroups = interactingAtoms.stream()
                        .map(Atom::getParentGroup)
                        .distinct()
                        .collect(Collectors.toList());

                interactingGroups.stream()
                        // drop interactions of a group with itself
                        .filter(interactingGroup -> !currentGroup.get().equals(interactingGroup))
                        .map(interactingGroup -> new PLIPInteraction(plipInteractionType, currentGroup.get(), interactingGroup))
                        .forEach(plipInteractions::add);
            }
        }

        return plipInteractions;
    }
}
