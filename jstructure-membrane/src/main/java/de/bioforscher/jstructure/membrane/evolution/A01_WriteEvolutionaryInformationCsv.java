package de.bioforscher.jstructure.membrane.evolution;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.feature.evolution.EvolutionaryInformation;
import de.bioforscher.jstructure.feature.interactions.PLIPInteractionContainer;
import de.bioforscher.jstructure.feature.interactions.PLIPIntraMolecularAnnotator;
import de.bioforscher.jstructure.feature.rigidity.BackboneRigidity;
import de.bioforscher.jstructure.feature.rigidity.DynaMineBridge;
import de.bioforscher.jstructure.feature.topology.OrientationsOfProteinsInMembranesAnnotator;
import de.bioforscher.jstructure.feature.topology.Topology;
import de.bioforscher.jstructure.membrane.MembraneConstants;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.jsoup.Jsoup;

import java.nio.file.Path;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

public class A01_WriteEvolutionaryInformationCsv {
    private static final Path DIRECTORY = MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY;
    private static final OrientationsOfProteinsInMembranesAnnotator ORIENTATIONS_OF_PROTEINS_IN_MEMBRANES_ANNOTATOR =
            new OrientationsOfProteinsInMembranesAnnotator();
    private static final DynaMineBridge DYNA_MINE_BRIDGE = new DynaMineBridge();
    private static final PLIPIntraMolecularAnnotator PLIP_INTRA_MOLECULAR_ANNOTATOR = new PLIPIntraMolecularAnnotator();

    public static void main(String[] args) {
        String header = "pdbId;chainId;resId;aa;topology;contact;total;rigidity;evol" + System.lineSeparator();

        String output = MembraneConstants.lines(DIRECTORY.resolve("ids.list"))
                .map(A01_WriteEvolutionaryInformationCsv::handleId)
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.joining(System.lineSeparator()));
        MembraneConstants.write(MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY.resolve("contact")
                        .resolve("aminoacids-evol-by-contact.csv"),
                header + output);
    }

    private static Optional<String> handleId(String id) {
        System.out.println(id);
        String pdbId = id.split("_")[0];
        String chainId = id.split("_")[1];

        try {
            Structure structure = StructureParser.source(DIRECTORY.resolve("pdb").resolve(pdbId + ".pdb"))
                    .minimalParsing(true)
                    .parse();
            Chain chain = structure.select()
                    .chainId(chainId)
                    .asChain();

            List<Double> evolutionaryScores = MembraneConstants.lines(DIRECTORY.resolve("evol").resolve(id + ".evol"))
                    .mapToDouble(Double::valueOf)
                    .boxed()
                    .collect(Collectors.toList());
            List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());
            if(aminoAcids.size() != evolutionaryScores.size()) {
                throw new IllegalArgumentException("sizes of lists do not match: " + aminoAcids.size() + " vs " + evolutionaryScores.size());
            }

            ORIENTATIONS_OF_PROTEINS_IN_MEMBRANES_ANNOTATOR.process(structure,
                    Jsoup.parse(MembraneConstants.lines(DIRECTORY.resolve("opm").resolve(pdbId + ".opm"))
                            .collect(Collectors.joining(System.lineSeparator()))));
            PLIP_INTRA_MOLECULAR_ANNOTATOR.process(chain,
                    Jsoup.parse(MembraneConstants.lines(DIRECTORY.resolve("plip").resolve(id + ".plip"))
                            .collect(Collectors.joining(System.lineSeparator()))));
            DYNA_MINE_BRIDGE.process(chain,
                    MembraneConstants.lines(DIRECTORY.resolve("dynamine").resolve(id + "_backbone.pred"))
                            .collect(Collectors.joining(System.lineSeparator())));

            for(int i = 0; i < aminoAcids.size(); i++) {
                aminoAcids.get(i).getFeatureContainer().addFeature(new EvolutionaryInformation(null, null, evolutionaryScores.get(i)));
            }

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
                    interactionContainer.getInteractions().size() + ";" +
                    StandardFormat.format(minMaxNormalize(aminoAcid.getFeature(BackboneRigidity.class).getBackboneRigidity(), min, max)) + ";" +
                    StandardFormat.format(aminoAcid.getFeature(EvolutionaryInformation.class).getInformation()));
        } catch (Exception e) {
            e.printStackTrace();
            return Optional.empty();
        }
    }

    private static double minMaxNormalize(double v, double min, double max) {
        return (v - min) / (max  - min);
    }
}
