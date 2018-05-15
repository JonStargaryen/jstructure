package de.bioforscher.jstructure.si.analysis;

import de.bioforscher.jstructure.efr.Start2FoldConstants;
import de.bioforscher.jstructure.efr.model.si.ResidueStructuralInformation;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.si.explorer.DataSource;
import de.bioforscher.jstructure.si.explorer.ExplorerChain;

import java.nio.file.Files;
import java.util.List;

/**
 * avg: spectrum b, red_white_green, minimum=-0.1, maximum=0.1
 * sum: spectrum b, red_white_green, minimum=-0.5, maximum=0.5
 * labels in inkscape: font: Fira Sans, bold, size: 32
 */
public class A01_WriteBFactorFiles {
    public static void main(String[] args) {
        DataSource.getInstance()
                .chains()
                .forEach(A01_WriteBFactorFiles::handleChain);
    }

    private static void handleChain(ExplorerChain explorerChain) {
        try {
            String entryId = explorerChain.getStfId();
            System.out.println(entryId);

            Chain chain = explorerChain.getChain();
            List<ResidueStructuralInformation> residueStructuralInformation = explorerChain.getResidues();

            // assign baseline: all atoms (including hetatms) have bfactor 0
            chain.atoms().forEach(atom -> atom.setBfactor(0));

            chain.aminoAcids()
                    .forEach(aminoAcid -> {
                        ResidueStructuralInformation information = residueStructuralInformation.get(aminoAcid.getAminoAcidIndex());
                        aminoAcid.atoms().forEach(atom -> {
                            atom.setBfactor((float) information.getAverageRmsdIncrease());
                        });
                    });
            Files.write(Start2FoldConstants.DATA_DIRECTORY.resolve("si").resolve("pymol").resolve(entryId + "-avg.pdb"),
                    chain.getPdbRepresentation().getBytes());

            chain.aminoAcids()
                    .forEach(aminoAcid -> {
                        ResidueStructuralInformation information = residueStructuralInformation.get(aminoAcid.getAminoAcidIndex());
                        aminoAcid.atoms().forEach(atom -> {
                            atom.setBfactor((float) information.getSumRmsdIncrease());
                        });
                    });
            Files.write(Start2FoldConstants.DATA_DIRECTORY.resolve("si").resolve("pymol").resolve(entryId + "-sum.pdb"),
                    chain.getPdbRepresentation().getBytes());
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
