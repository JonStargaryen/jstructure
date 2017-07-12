package de.bioforscher.jstructure.mutation.impl;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Element;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.ProteinParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.junit.Assert;
import org.junit.Test;

/**
 * Test for the residue mutator.
 * Created by bittrich on 7/6/17.
 */
public class ResidueMutatorServiceImplTest {
    @Test
    public void shouldIntroduceMutation() {
        String chainId = "A";
        int position = 27;
        AminoAcid.Family target = AminoAcid.Family.ALANINE;

        Protein protein = ProteinParser.source("2lzm").parse();
        Protein mutatedProtein = new ResidueMutatorServiceImpl().mutateResidue(protein, chainId, position, target);

        AminoAcid mutatedGroup = mutatedProtein.select()
                .chainName(chainId)
                .residueNumber(position)
                .asAminoAcid();
        AminoAcid previousGroup = mutatedGroup.getPreviousAminoAcid().get();
        AminoAcid nextGroup = mutatedGroup.getNextAminoAcid().get();

        // check for correct name
        Assert.assertEquals("group name does not match", target.getGroupPrototype().getThreeLetterCode(), mutatedGroup.getThreeLetterCode());
        // check for correct atom count - subtract 1 as prototype atoms contains OXT
        Assert.assertEquals("atom count does not match", target.getGroupPrototype()
                .getPrototypeAtoms()
                .stream()
                .map(Atom::getElement)
                .filter(Element::isHeavyAtom)
                .count() - 1,
                mutatedGroup.getAtoms().size());
        // check for spatial alignment
        Assert.assertTrue("groups not spatially neighbored", 4.0 < mutatedGroup.calculate()
                .center()
                .distanceFast(previousGroup.calculate().center()));
        Assert.assertTrue("groups not spatially neighbored", 4.0 < mutatedGroup.calculate()
                .center()
                .distanceFast(nextGroup.calculate().center()));
    }
}