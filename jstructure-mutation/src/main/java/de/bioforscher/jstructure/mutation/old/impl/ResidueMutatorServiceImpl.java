package de.bioforscher.jstructure.mutation.old.impl;

import de.bioforscher.jstructure.align.AlignmentPolicy;
import de.bioforscher.jstructure.align.StructureAlignmentBuilder;
import de.bioforscher.jstructure.align.impl.SingleValueDecompositionAligner;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.identifier.ChainIdentifier;
import de.bioforscher.jstructure.model.structure.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.mutation.old.ResidueMutatorService;

import java.util.List;

/**
 * Impl of the residue mutator.
 * Created by bittrich on 7/6/17.
 */
@Deprecated
public class ResidueMutatorServiceImpl implements ResidueMutatorService {
    private final SingleValueDecompositionAligner singleValueDecompositionAligner;

    public ResidueMutatorServiceImpl() {
        this.singleValueDecompositionAligner = new SingleValueDecompositionAligner();
    }

    @Override
    public Protein mutateResidue(Protein originalProtein, String chainId, int residueNumber, AminoAcid.Family targetAminoAcid) {
        try {
            GroupContainer originalResidueContainer = originalProtein.select()
                    .chainName(chainId)
                    .residueNumber(residueNumber)
                    .asGroupContainer();
            Group originalGroup = originalResidueContainer.getGroups().get(0);

            // deep-copy original protein
            Protein mutatedProtein = new Protein(originalProtein);
            // select chain to mutate
            Chain mutatedChain = mutatedProtein.select()
                    .chainName(chainId)
                    .asChain();
            Group positionInMutatedChain = mutatedChain.select()
                    .residueNumber(residueNumber)
                    .asGroup();

            // get prototype atoms from target amino acid - we need GroupContainer instances for the superimposition
            Group mutatedResidue = AminoAcid.Family.createAminoAcid(targetAminoAcid.getThreeLetterCode(), originalGroup.getResidueIdentifier(), originalGroup.isLigand());

            // assign prototype atoms
            targetAminoAcid.getGroupPrototype()
                    .getPrototypeAtoms()
                    .stream()
                    // ignore hydrogen atoms and terminal oxygen
                    .filter(atom -> atom.getElement().isHeavyAtom() && !atom.getName().equals("OXT"))
                    // clone atoms, so that original atoms are by no means compromised
                    .map(Atom::new)
                    .forEach(mutatedResidue::addAtom);

            Chain mutatedResidueContainer = new Chain(ChainIdentifier.UNKNOWN_CHAIN_ID);
            mutatedResidueContainer.addGroup(mutatedResidue);
            // transform coordinates
            StructureAlignmentBuilder.StructureAlignmentStep query = StructureAlignmentBuilder.builder(originalResidueContainer, mutatedResidueContainer)
                    .matchingBehavior(AlignmentPolicy.MatchingBehavior.comparableAtomNames)
                    .manipulationBehavior(AlignmentPolicy.ManipulationBehavior.COPY);
            singleValueDecompositionAligner.align(query)
                    .getTransformation()
                    .transform(mutatedResidue);

            // replace old residue with mutated one
            mutatedChain.getGroups().remove(positionInMutatedChain);
            mutatedChain.addGroup(mutatedResidue);

            // minimize energy or select reasonable rotamer
            //TODO optionally impl

            // renumber, so that pdbSerials are consistent
            renumber(mutatedProtein);

            return mutatedProtein;
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    private void renumber(Protein protein) {
        int pdbSerial = 1;
        for(Chain chain : protein.getChains()) {
            List<Group> groups = chain.getGroups();
            // sort groups by residue numbers
            groups.sort(ResidueIdentifier.GROUP_COMPARATOR);
            for(Group group : groups) {
                for(Atom atom : group.getAtoms()) {
                    atom.setPdbSerial(pdbSerial);
                    pdbSerial++;
                }
            }
        }
    }
}
