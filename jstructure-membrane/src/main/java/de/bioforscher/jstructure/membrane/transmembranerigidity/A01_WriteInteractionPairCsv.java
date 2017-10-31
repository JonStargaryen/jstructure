package de.bioforscher.jstructure.membrane.transmembranerigidity;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.feature.rigidity.BackboneRigidity;
import de.bioforscher.jstructure.feature.rigidity.DynaMineBridge;
import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.feature.sse.dssp.DictionaryOfProteinSecondaryStructure;
import de.bioforscher.jstructure.feature.topology.OrientationsOfProteinsInMembranesAnnotator;
import de.bioforscher.jstructure.feature.topology.Topology;
import de.bioforscher.jstructure.membrane.MembraneConstants;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.jsoup.Jsoup;

import java.nio.file.Path;
import java.util.Optional;
import java.util.stream.Collectors;

public class A01_WriteInteractionPairCsv {
    private static final DictionaryOfProteinSecondaryStructure DICTIONARY_OF_PROTEIN_SECONDARY_STRUCTURE = new DictionaryOfProteinSecondaryStructure();
    private static final DynaMineBridge DYNA_MINE_BRIDGE = new DynaMineBridge();
    private static final OrientationsOfProteinsInMembranesAnnotator ORIENTATIONS_OF_PROTEINS_IN_MEMBRANES_ANNOTATOR = new OrientationsOfProteinsInMembranesAnnotator();

    public static void main(String[] args) {
        String header = "pdbId;resNum;sse;sseSize;topology;rigidity" + System.lineSeparator();

        String alpha = composeCsv(MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY);
        MembraneConstants.write(MembraneConstants.DATASETS_DIRECTORY.resolve("transmembranerigidity").resolve("alpha.csv"), header + alpha);

        String beta = composeCsv(MembraneConstants.PDBTM_NR_BETA_DATASET_DIRECTORY);
        MembraneConstants.write(MembraneConstants.DATASETS_DIRECTORY.resolve("transmembranerigidity").resolve("beta.csv"), header + beta);
    }

    private static String composeCsv(Path directory) {
        return MembraneConstants.lines(directory.resolve("ids.list"))
                .map(line -> handleLine(line, directory))
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.joining(System.lineSeparator()));
    }

    private static Optional<String> handleLine(String line, Path directory) {
        try {
            System.out.println(line);
            String pdbId = line.split("_")[0];
            String chainId = line.split("_")[1];

            Structure structure = StructureParser.source(directory.resolve("pdb").resolve(pdbId + ".pdb"))
                    .minimalParsing(true)
                    .parse();
            Chain chain = structure.select()
                    .chainName(chainId)
                    .asChain();

            DICTIONARY_OF_PROTEIN_SECONDARY_STRUCTURE.process(structure);
            DYNA_MINE_BRIDGE.process(chain, MembraneConstants.lines(directory.resolve("dynamine").resolve(line + "_backbone.pred"))
                    .collect(Collectors.joining(System.lineSeparator())));
            ORIENTATIONS_OF_PROTEINS_IN_MEMBRANES_ANNOTATOR.process(structure, Jsoup.parse(MembraneConstants.lines(directory.resolve("opm").resolve(pdbId + ".opm"))
                    .collect(Collectors.joining(System.lineSeparator()))));

            double min = chain.aminoAcids()
                    .mapToDouble(aa -> aa.getFeature(BackboneRigidity.class).getBackboneRigidity())
                    .min()
                    .getAsDouble();
            double max = chain.aminoAcids()
                    .mapToDouble(aa -> aa.getFeature(BackboneRigidity.class).getBackboneRigidity())
                    .max()
                    .getAsDouble();

            return Optional.of(chain.aminoAcids()
                    .map(aminoAcid -> composeLine(aminoAcid, min, max))
                    .collect(Collectors.joining(System.lineSeparator())));
        } catch (Exception e) {
            e.printStackTrace();
            return Optional.empty();
        }
    }

    private static String composeLine(AminoAcid aminoAcid, double min, double max) {
        GenericSecondaryStructure sse = aminoAcid.getFeature(GenericSecondaryStructure.class);
        String aaSse = sse.getSecondaryStructure().getReducedRepresentation();
        int sizeOfSse = sse.getSurroundingSecondaryStructureElement(aminoAcid).getSize();

        return aminoAcid.getParentChain().getChainIdentifier().getFullName() + ";" +
                aminoAcid.getResidueIdentifier() + ";" +
                aaSse + ";" +
                sizeOfSse + ";" +
                (aminoAcid.getFeature(Topology.class).isTransmembrane() ? "I" : "o") + ";" +
                StandardFormat.format(minMaxNormalize(aminoAcid.getFeature(BackboneRigidity.class).getBackboneRigidity(), min, max));
    }

    private static double minMaxNormalize(double v, double min, double max) {
        return (v - min) / (max  - min);
    }
}
