package de.bioforscher.jstructure.model.structure.container;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Protein;

import java.util.List;
import java.util.stream.Stream;

/**
 * Specifies the capabilities of a getChain container (mostly a {@link Protein}).
 * Created by S on 30.09.2016.
 */
public interface ChainContainer extends GroupContainer {
    /**
     * Never manipulate the returned collection as it is not guaranteed the actually modify the internal list(s).
     * @return all associated chains
     */
    List<Chain> getChains();

    /**
     * Access to all chain objects associated to this container.
     * @return a stream of chains
     */
    default Stream<Chain> chains() {
        return getChains().stream();
    }

    /**
     * Access to all chains which actually contains amino acids.
     * @return a stream of chains with at least 1 amino acid
     */
    default Stream<Chain> chainsWithAminoAcids() {
        return chains()
                .filter(chain -> chain.aminoAcids().count() > 0);
    }

    //TODO do Protein instances really return groups/atoms/sequence for chains beyond the first?
//    @Override
//    default List<Group> getGroups() {
//        return chains()
//                .flatMap(Chain::groups)
//                .collect(Collectors.toList());
//    }
//
//    @Override
//    default List<Atom> getAtoms() {
//        return chains()
//                .flatMap(Chain::groups)
//                .flatMap(Group::atoms)
//                .collect(Collectors.toList());
//    }
//
//    @Override
//    default Stream<AminoAcid> aminoAcids() {
//        return chains()
//                .flatMap(Chain::groups)
//                .filter(Group::isAminoAcid)
//                .map(AminoAcid.class::cast);
//    }
//
//    @Override
//    default String getAminoAcidSequence() {
//        return aminoAcids()
//                .map(AminoAcid::getOneLetterCode)
//                .collect(Collectors.joining());
//    }

    @Override
    default ChainContainer createCopy() {
        return (ChainContainer) GroupContainer.super.createCopy();
    }
}
