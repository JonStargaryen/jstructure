package de.bioforscher.jstructure.si.analysis;

import de.bioforscher.jstructure.efr.Start2FoldConstants;
import de.bioforscher.jstructure.efr.model.si.ContactStructuralInformation;
import de.bioforscher.jstructure.mathematics.SetOperations;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.si.explorer.DataSource;
import de.bioforscher.jstructure.si.explorer.ExplorerChain;

import java.nio.file.Files;
import java.util.List;
import java.util.Optional;
import java.util.StringJoiner;
import java.util.stream.Collectors;

public class A02_WriteHeatmapFiles {
    public static void main(String[] args) {
        DataSource.getInstance()
                .chains()
                .forEach(A02_WriteHeatmapFiles::handleChain);
    }

    private static void handleChain(ExplorerChain explorerChain) {
        try {
            String entryId = explorerChain.getStfId();
            System.out.println(entryId);

            Chain chain = explorerChain.getChain();
            List<ContactStructuralInformation> contactStructuralInformation = explorerChain.getContacts();

            StringJoiner stringJoiner = new StringJoiner(System.lineSeparator(),
                    "res1,res2,avg" + System.lineSeparator(),
                    "");
            List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());
            SetOperations.cartesianProductOf(aminoAcids, aminoAcids)
                    .forEach(pair -> {
                        int residueIdentifier1 = pair.getLeft().getResidueIdentifier().getResidueNumber();
                        int residueIdentifier2 = pair.getRight().getResidueIdentifier().getResidueNumber();
                        Optional<ContactStructuralInformation> information = contactStructuralInformation.stream()
                                .filter(contact -> (contact.getResidueIdentifier1() == residueIdentifier1 && contact.getResidueIdentifier2() == residueIdentifier2) ||
                                        (contact.getResidueIdentifier1() == residueIdentifier2 && contact.getResidueIdentifier2() == residueIdentifier1))
                                .findFirst();
                        if(information.isPresent()) {
                            stringJoiner.add(residueIdentifier1 + "," +
                                    residueIdentifier2 + "," +
                                    information.get().getAverageRmsdIncrease());
                        } else {
                            stringJoiner.add(residueIdentifier1 + "," +
                                    residueIdentifier2 + "," +
                                    "0.0");
                        }
                    });
            Files.write(Start2FoldConstants.DATA_DIRECTORY.resolve("si").resolve("heatmaps").resolve(entryId + ".csv"),
                    stringJoiner.toString().getBytes());
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
