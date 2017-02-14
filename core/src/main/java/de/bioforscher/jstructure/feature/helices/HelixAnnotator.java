package de.bioforscher.jstructure.feature.helices;

import de.bioforscher.jstructure.feature.sse.SecStrucState;
import de.bioforscher.jstructure.feature.sse.SecondaryStructureAnnotator;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.selection.IntegerRange;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.List;

/**
 * Annotate properties of helices.
 * Created by bittrich on 2/13/17.
 */
@FeatureProvider(provides = HelixAnnotator.HELIX_PROPERTIES, requires = SecondaryStructureAnnotator.SECONDARY_STRUCTURE_STATES)
public class HelixAnnotator extends AbstractFeatureProvider {
    private static final Logger logger = LoggerFactory.getLogger(HelixAnnotator.class);
    public static final String HELIX_PROPERTIES = "HELIX_PROPERTIES";

    @Override
    protected void processInternally(Protein protein) {
        List<Helix> globalListOfHelices = new ArrayList<>();
        for(Chain chain : protein.getChains()) {
            Group startGroupOfHelix = null;
            Group lastGroupOfHelix = null;
            for(Group group : chain.getGroups()) {
                if (!group.isAminoAcid()) {
                    continue;
                }

                boolean isHelix = isHelix(group);
                // do nothing if this residue is not in a helix and it does not cause a helix to stop
                if (!isHelix && lastGroupOfHelix == null) {
                    continue;
                }

                // helix keeps extending or just started
                if (isHelix) {
                    if (startGroupOfHelix != null) {
                        lastGroupOfHelix = group;
                    } else {
                        startGroupOfHelix = group;
                    }
                    continue;
                }

                // allow 1 and 2 non-helix residue and still combine fragments into 1 helix
                // 'Any helical segments separated by only one or two residues were combined.'
                // Law, 2016 - Examining the conservation of kinks in alpha helices
                if(group.getResidueNumber() - lastGroupOfHelix.getResidueNumber() <= 2) {
                    continue;
                }

                // helix just ended
                Helix helix = new Helix(startGroupOfHelix, lastGroupOfHelix);
                Selection.on(chain)
                        .residueNumber(new IntegerRange(startGroupOfHelix.getResidueNumber(), lastGroupOfHelix.getResidueNumber()))
                        .asFilteredGroups()
                        .forEach(groupToAssignTo ->
                            groupToAssignTo.setFeature(HELIX_PROPERTIES, helix)
                        );
                logger.debug("annotated helix: {}", helix);
            }
        }

        protein.setFeature(HELIX_PROPERTIES, globalListOfHelices);
    }

    private boolean isHelix(Group group) {
        return group.getFeatureAsList(SecStrucState.class, SecondaryStructureAnnotator.SECONDARY_STRUCTURE_STATES).stream()
                .anyMatch(secStrucState -> secStrucState.getSecondaryStructure().isHelixType());
    }
}
