package de.bioforscher.jstructure.mutation.impl;

import de.bioforscher.jstructure.align.AlignmentPolicy;
import de.bioforscher.jstructure.align.StructureAlignmentQuery;
import de.bioforscher.jstructure.align.impl.SingleValueDecompositionAligner;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.identifier.ChainIdentifier;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.mutation.MutatorService;

import java.util.List;

/**
 * Impl of the simple residue mutator without any refinement of the side-chain location.
 * Created by bittrich on 7/13/17.
 */
public class SimpleMutatorServiceImpl implements MutatorService {
    private final SingleValueDecompositionAligner singleValueDecompositionAligner;

    public SimpleMutatorServiceImpl() {
        this.singleValueDecompositionAligner = new SingleValueDecompositionAligner();
    }

    @Override
    public Structure mutateAminoAcid(Structure originalProtein,
                                     ChainIdentifier chainIdentifier,
                                     ResidueIdentifier aminoAcidToMutate,
                                     AminoAcid.Family targetAminoAcid) {
        try {
            GroupContainer originalResidueContainer = originalProtein.select()
                    .chainName(chainIdentifier.getChainId())
                    .residueNumber(aminoAcidToMutate.getResidueNumber())
                    .asIsolatedStructure();
            Group originalGroup = originalResidueContainer.getGroups().get(0);

            // deep-copy original protein
            Structure mutatedProtein = originalProtein.createDeepCopy();
            // select chain to mutate
            Chain mutatedChain = mutatedProtein.select()
                    .chainName(chainIdentifier.getChainId())
                    .asChain();
            Group positionInMutatedChain = mutatedChain.select()
                    .residueNumber(aminoAcidToMutate.getResidueNumber())
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
                    .map(Atom::createDeepCopy)
                    .forEach(mutatedResidue::addAtom);

            Chain mutatedResidueContainer = new Chain(ChainIdentifier.UNKNOWN_CHAIN_IDENTIFIER);
            mutatedResidueContainer.addGroup(mutatedResidue);
            // transform coordinates
            StructureAlignmentQuery query = StructureAlignmentQuery.of(originalResidueContainer, mutatedResidueContainer)
                    .matchingBehavior(AlignmentPolicy.MatchingBehavior.comparableAtomNames);
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

    private void renumber(Structure protein) {
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