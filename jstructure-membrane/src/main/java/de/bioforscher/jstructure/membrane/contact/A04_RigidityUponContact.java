package de.bioforscher.jstructure.membrane.contact;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.interactions.PLIPIntraMolecularAnnotator;
import de.bioforscher.jstructure.feature.rigidity.BackboneRigidity;
import de.bioforscher.jstructure.feature.rigidity.DynaMineBridge;
import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
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

public class A04_RigidityUponContact {
    private static final OrientationsOfProteinsInMembranesAnnotator ORIENTATIONS_OF_PROTEINS_IN_MEMBRANES_ANNOTATOR =
            new OrientationsOfProteinsInMembranesAnnotator();
    private static final PLIPIntraMolecularAnnotator PLIP_INTRA_MOLECULAR_ANNOTATOR =
            new PLIPIntraMolecularAnnotator();
    private static final DynaMineBridge DYNA_MINE_BRIDGE = new DynaMineBridge();
    private static final Path DIRECTORY = MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY;

    public static void main(String[] args) {
        String header = "pdbId;chainId;resId;aa;topology;contact;rigidity" + System.lineSeparator();

        String output = MembraneConstants.lines(DIRECTORY.resolve("ids.list"))
                .map(A04_RigidityUponContact::handleLine)
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.joining(System.lineSeparator()));
        MembraneConstants.write(MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY.resolve("contact")
                        .resolve("aminoacids-rigidity-by-contact.csv"),
                header + output);
    }

    private static Optional<String> handleLine(String line) {
        try {
            System.out.println(line);
            String pdbId = line.split("_")[0];
            String chainId = line.split("_")[1];

            Structure structure = StructureParser.source(DIRECTORY.resolve("pdb").resolve(pdbId + ".pdb"))
                    .minimalParsing(true)
                    .parse();
            Chain chain = structure.select()
                    .chainName(chainId)
                    .asChain();

            ORIENTATIONS_OF_PROTEINS_IN_MEMBRANES_ANNOTATOR.process(structure,
                    Jsoup.parse(MembraneConstants.lines(DIRECTORY.resolve("opm").resolve(pdbId + ".opm"))
                            .collect(Collectors.joining(System.lineSeparator()))));
            PLIP_INTRA_MOLECULAR_ANNOTATOR.process(chain,
                    Jsoup.parse(MembraneConstants.lines(DIRECTORY.resolve("plip").resolve(line + ".plip"))
                            .collect(Collectors.joining(System.lineSeparator()))));
            DYNA_MINE_BRIDGE.process(chain,
                    MembraneConstants.lines(DIRECTORY.resolve("dynamine").resolve(line + "_backbone.pred"))
                            .collect(Collectors.joining(System.lineSeparator())));

            double min = chain.aminoAcids()
                    .mapToDouble(aa -> aa.getFeature(BackboneRigidity.class).getBackboneRigidity())
                    .min()
                    .getAsDouble();
            double max = chain.aminoAcids()
                    .mapToDouble(aa -> aa.getFeature(BackboneRigidity.class).getBackboneRigidity())
                    .max()
                    .getAsDouble();

            String output = chain.aminoAcids()
                    .map(aminoAcid -> handleAminoAcid(pdbId, chainId, aminoAcid, min, max))
                    .filter(Optional::isPresent)
                    .map(Optional::get)
                    .collect(Collectors.joining(System.lineSeparator()));
            return Optional.of(output);
        } catch (Exception e) {
            e.printStackTrace();
            return Optional.empty();
        }
    }

    private static Optional<String> handleAminoAcid(String pdbId,
                                                    String chainId,
                                                    AminoAcid aminoAcid,
                                                    double min,
                                                    double max) {
        try {
            boolean tm = aminoAcid.getFeature(Topology.class).isTransmembrane();
            GenericSecondaryStructure.SecondaryStructureElement surroundingSecondaryStructureElement =
                    aminoAcid.getFeature(GenericSecondaryStructure.class).getSurroundingSecondaryStructureElement(aminoAcid);
            // skip non-transmembrane parts when there is no helix or it is too short
            if (!tm && (!surroundingSecondaryStructureElement.getReducedType().equals("H") ||
                    surroundingSecondaryStructureElement.getSize() <= MembraneConstants.MINIMUM_HELIX_LENGTH)) {
                return Optional.empty();
            }

            // evaluate interaction types
            PLIPInteractionContainer interactionContainer = aminoAcid.getFeature(PLIPInteractionContainer.class);
            boolean inContact = tm ?
                    interactionContainer.getInteractions()
                            .stream()
                            .anyMatch(MembraneConstants::isTmHelixInteraction) :
                    interactionContainer.getInteractions()
                            .stream()
                            .anyMatch(MembraneConstants::isNtmHelixInteraction);

            return Optional.of(pdbId + ";" +
                    chainId + ";" +
                    aminoAcid.getResidueIdentifier() + ";" +
                    aminoAcid.getOneLetterCode() + ";" +
                    (tm ? "I" : "o") + ";" +
                    (inContact ? "contact" : "no-contact") + ";" +
                    StandardFormat.format(minMaxNormalize(aminoAcid.getFeature(BackboneRigidity.class).getBackboneRigidity(), min, max)));
        } catch (Exception e) {
            e.printStackTrace();
            return Optional.empty();
        }
    }

    private static double minMaxNormalize(double v, double min, double max) {
        return (v - min) / (max  - min);
    }
}
