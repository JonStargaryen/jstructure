package studies.plipvalidation;

import de.bioforscher.jstructure.feature.interactions.PLIPAnnotator;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.sse.dssp.DSSPSecondaryStructure;
import de.bioforscher.jstructure.feature.sse.dssp.DictionaryOfProteinSecondaryStructure;
import de.bioforscher.jstructure.model.feature.FeatureContainer;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.ProteinParser;

/**
 * Created by bittrich on 6/29/17.
 */
public class PLIPValidation {
    public static void main(String[] args) {
        process("1acj");
    }

    private static void process(String pdbId) {
        Protein protein = ProteinParser.source(pdbId).parse();

        try {
            new DictionaryOfProteinSecondaryStructure().process(protein);
            new PLIPAnnotator().process(protein);

            for (Chain chain : protein.getChains()) {
                for (Group group : chain.getGroups()) {
                    if (!(group instanceof AminoAcid)) {
                        continue;
                    }

                    FeatureContainer featureContainer = group.getFeatureContainer();
                    PLIPInteractionContainer plipInteractionContainer = featureContainer.getFeature(PLIPInteractionContainer.class);
                    DSSPSecondaryStructure dsspSecondaryStructure = featureContainer.getFeature(DSSPSecondaryStructure.class);

                    System.out.println(group);
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}