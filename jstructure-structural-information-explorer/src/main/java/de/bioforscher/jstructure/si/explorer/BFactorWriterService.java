package de.bioforscher.jstructure.si.explorer;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.si.explorer.model.ResidueStructuralInformation;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.stream.Collectors;

public class BFactorWriterService {
    private static final BFactorWriterService INSTANCE = new BFactorWriterService();

    private BFactorWriterService() {

    }

    public static BFactorWriterService getInstance() {
        return INSTANCE;
    }

    public void writePDBFileWithBFactors(Chain chain,
                                         List<ResidueStructuralInformation> structuralInformationList,
                                         Path outputPathAverage,
                                         Path outputPathMax) {
        try {
            List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());
            for (int i = 0; i < chain.aminoAcids().count(); i++) {
                AminoAcid aminoAcid = aminoAcids.get(i);
                ResidueStructuralInformation structuralInformation = structuralInformationList.get(i);
                aminoAcid.atoms().forEach(atom -> atom.setBfactor((float) structuralInformation.getAverageRmsdIncrease()));
            }
            Files.write(outputPathAverage, chain.getPdbRepresentation().getBytes());

            for (int i = 0; i < chain.aminoAcids().count(); i++) {
                AminoAcid aminoAcid = aminoAcids.get(i);
                ResidueStructuralInformation structuralInformation = structuralInformationList.get(i);
                aminoAcid.atoms().forEach(atom -> atom.setBfactor((float) structuralInformation.getMaximumRmsdIncrease()));
            }
            Files.write(outputPathMax, chain.getPdbRepresentation().getBytes());
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }
}
