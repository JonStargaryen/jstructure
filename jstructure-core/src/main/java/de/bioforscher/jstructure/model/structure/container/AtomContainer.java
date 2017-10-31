package de.bioforscher.jstructure.model.structure.container;

import de.bioforscher.jstructure.mathematics.LinearAlgebra;
import de.bioforscher.jstructure.mathematics.Calculable;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.selection.Selectable;

import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Defines the capability to convert itself to <tt>ATOM</tt> records according to <tt>PDB</tt> format.
 * Created by S on 28.09.2016.
 */
public interface AtomContainer extends StructureContainer, Selectable,
        Calculable<LinearAlgebra.AtomContainerLinearAlgebra> {
    /**
     * Never manipulate the returned collection as it is not guaranteed to actually modify the internal list(s).
     * @return all associated atoms
     */
    List<Atom> getAtoms();

    /**
     * Access to all atoms of this container.
     * @return a select of all atoms in this container
     */
    default Stream<Atom> atoms() {
        return getAtoms().stream();
    }

    default String getPdbRepresentation() {
        Structure protein;
        try {
            protein = getAtoms().get(0).getParentGroup().getParentChain().getParentStructure();
        } catch (Exception e) {
            protein = Structure.UNKNOWN_STRUCTURE;
        }

        // ensure ordering of atoms according to their pdbSerials
        List<Atom> sortedAtoms = atoms()
                .sorted(Comparator.comparingInt(Atom::getPdbSerial))
                .collect(Collectors.toList());

        return getPdbRepresentation(sortedAtoms, protein);
    }

    static String getPdbRepresentation(List<Atom> sortedAtoms, Structure protein) {
        StringBuilder stringBuilder = new StringBuilder();
        stringBuilder.append(protein.getHeader());

        int previousPdbSerial = 0;
        Group previousGroup = null;
        Chain previousChain = null;
        Set<Chain> terminatedChains = new HashSet<>();

        for(Atom atom : sortedAtoms) {
            Group currentGroup = atom.getParentGroup();
            Chain currentChain = currentGroup.getParentChain();
            boolean currentGroupIsInTerminatedChain = currentGroup.isLigand();
            boolean previousChainIsNotMarkedAsTerminated = !terminatedChains.contains(previousChain);
            boolean currentChainIsNotMarkedAsTerminated = !terminatedChains.contains(currentChain);

            // chain changed and previous chain was not yet terminated
            if(previousChain != null) {
                if (!currentChain.equals(previousChain) && previousChainIsNotMarkedAsTerminated) {
                    terminatedChains.add(previousChain);
                    writeTerminateRecord(stringBuilder, previousChain, previousGroup, previousPdbSerial);
                    // the current group is a ligand, but the current chain was not yet terminated
                } else if (currentGroupIsInTerminatedChain && currentChainIsNotMarkedAsTerminated) {
                    // check if chain needs TER record (only polymer chains require it)
                    if(currentChain.groups().count() != currentChain.ligands().count()) {
                        terminatedChains.add(currentChain);
                        writeTerminateRecord(stringBuilder, currentChain, previousGroup, previousPdbSerial);
                    }
                }
            }
            previousPdbSerial = atom.getPdbSerial();
            previousGroup = currentGroup;
            previousChain = currentChain;

            // print ATOM record as usual
            stringBuilder.append(atom.getPdbRepresentation())
                    .append(System.lineSeparator());
        }

        stringBuilder.append("END")
                .append(System.lineSeparator());

        return stringBuilder.toString();
    }

    static void writeTerminateRecord(StringBuilder stringBuilder, Chain chain, Group previousGroup, int previousPdbSerial) {
        // TER     961      ASP A  62
        stringBuilder.append("TER   ")
                .append(String.format("%5d", (previousPdbSerial + 1)))
                .append("      ")
                .append(previousGroup.getThreeLetterCode())
                .append(" ")
                .append(chain.getChainIdentifier().getChainId())
                .append(String.format("%4d", previousGroup.getResidueIdentifier().getResidueNumber()))
                .append(previousGroup.getResidueIdentifier().getInsertionCode())
                .append(System.lineSeparator());
    }
}