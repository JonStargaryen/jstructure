package de.bioforscher.jstructure.mutation.impl;

import de.bioforscher.jstructure.align.AlignmentPolicy;
import de.bioforscher.jstructure.align.StructureAligner;
import de.bioforscher.jstructure.model.structure.*;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.identifier.ChainIdentifier;
import de.bioforscher.jstructure.model.structure.identifier.ResidueIdentifier;

import java.util.Comparator;
import java.util.List;

/**
 * Introduces a mutation for a given residue in a protein structure.
 * Created by bittrich on 7/6/17.
 */
class ResidueMutatorService {
    /**
     * Mutate a given position in a protein.
     * @param originalProtein the original protein which will not be manipulated
     * @param chainId the chainId of the mutation
     * @param residueNumber the residueNumber of the mutation
     * @param targetAminoAcid the desired target amino acid
     * @return a new {@link Protein} instance featuring the mutation
     */
    public Protein mutate(Protein originalProtein, String chainId, int residueNumber, AminoAcid.Family targetAminoAcid) {
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
            Group mutatedResidue = targetAminoAcid.getRepresentingClass().getConstructor(ResidueIdentifier.class, boolean.class)
                    .newInstance(originalGroup.getResidueIdentifier(), originalGroup.isLigand());

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
            StructureAligner.builder(originalResidueContainer, mutatedResidueContainer)
                    .matchingBehavior(AlignmentPolicy.MatchingBehavior.comparableAtomNames)
                    .manipulationBehavior(AlignmentPolicy.ManipulationBehavior.COPY)
                    .align()
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
            //TODO maybe provide standard comparator via ResidueIdentifier class
            groups.sort(Comparator.comparingInt(group -> group.getResidueIdentifier().getResidueNumber()));
            for(Group group : groups) {
                for(Atom atom : group.getAtoms()) {
                    atom.setPdbSerial(pdbSerial);
                    pdbSerial++;
                }
            }
        }
    }
}
