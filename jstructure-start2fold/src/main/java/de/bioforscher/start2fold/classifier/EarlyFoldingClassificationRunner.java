package de.bioforscher.start2fold.classifier;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.feature.asa.AccessibleSurfaceAreaCalculator;
import de.bioforscher.jstructure.feature.energyprofile.EgorAgreementCalculator;
import de.bioforscher.jstructure.feature.interactions.PLIPIntraMolecularAnnotator;
import de.bioforscher.jstructure.feature.interactions.PLIPRestServiceQuery;
import de.bioforscher.jstructure.feature.loopfraction.LoopFraction;
import de.bioforscher.jstructure.feature.loopfraction.LoopFractionCalculator;
import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.model.identifier.ChainIdentifier;
import de.bioforscher.jstructure.model.identifier.IdentifierFactory;
import de.bioforscher.jstructure.model.identifier.ProteinIdentifier;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.aminoacid.Proline;
import de.bioforscher.start2fold.model.vector.RawFeatureVector;
import de.bioforscher.start2fold.model.vector.SmoothedFeatureVector;
import org.jsoup.nodes.Document;
import weka.classifiers.Classifier;
import weka.core.*;

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

public class EarlyFoldingClassificationRunner {
    private final EgorAgreementCalculator EGOR_AGREEMENT_CALCULATOR = new EgorAgreementCalculator();
    private final LoopFractionCalculator LOOP_FRACTION_CALCULATOR = new LoopFractionCalculator();
    private final AccessibleSurfaceAreaCalculator ACCESSIBLE_SURFACE_AREA_CALCULATOR = new AccessibleSurfaceAreaCalculator();
    private final PLIPIntraMolecularAnnotator PLIP_INTRA_MOLECULAR_ANNOTATOR = new PLIPIntraMolecularAnnotator();
    private final Classifier model;
    private static final EarlyFoldingClassificationRunner INSTANCE = new EarlyFoldingClassificationRunner();

    private EarlyFoldingClassificationRunner() {
        try {
            InputStream inputStream = Thread.currentThread()
                    .getContextClassLoader()
                    .getResourceAsStream("model/efr-classifier.model");
            model = (Classifier) SerializationHelper.read(inputStream);
        } catch (Exception e) {
            System.out.println("could not load GMLVQ-model");
            e.printStackTrace();
            throw new RuntimeException();
        }
    }

    public static EarlyFoldingClassificationRunner getInstance() {
        return INSTANCE;
    }

    public void process(Structure structure, Path outputPath) throws IOException {
        // report structure characteristics
        System.out.println("structure: " + structure.getProteinIdentifier().getFullName() + "\n" +
                "chains: " + structure.chainsWithAminoAcids()
                .map(Chain::getChainIdentifier)
                .map(ChainIdentifier::getChainId)
                .collect(Collectors.toList()) + "\n" +
                "total residues: " + structure.aminoAcids().count());

        System.out.println();

        // compute features
        System.out.println("computing residue-level features");

        // start with PLIP to fail fast
        System.out.println("querying PLIP-REST-Service");
        try {
            // try to annotate by standard routine
            PLIP_INTRA_MOLECULAR_ANNOTATOR.process(structure);
            System.out.println("fetched PLIP contacts");
        } catch (Exception e1) {
            try {
                // potential non-pdb-entry, try to compute on-the-fly
                structure.chainsWithAminoAcids().forEach(chain -> {
                    Document document = PLIPRestServiceQuery.calculateIntraChainDocument(chain);
                    PLIP_INTRA_MOLECULAR_ANNOTATOR.process(chain, document);
                });
                System.out.println("computed PLIP contacts");
            } catch (Exception e2) {
                System.out.println("failed: could not compute PLIP contacts");
                e2.printStackTrace();
                return;
            }
        }

        System.out.println("computing energy profiles");
        EGOR_AGREEMENT_CALCULATOR.process(structure);

        System.out.println("annotating secondary structure elements");
        LOOP_FRACTION_CALCULATOR.process(structure);

        System.out.println("computing relative accessible surface area");
        ACCESSIBLE_SURFACE_AREA_CALCULATOR.process(structure);

        // assign feature vectors
        structure.aminoAcids().forEach(RawFeatureVector::assignRawFeatureVector);

        // smooth feature vectors
        structure.chainsWithAminoAcids().forEach(chain -> {
            List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());
            aminoAcids.forEach(aminoAcid -> {
                SmoothedFeatureVector.assignSmoothedFeatureVector(aminoAcids, aminoAcid);
            });
        });

        // classify each residue
        StringJoiner outputJoiner = new StringJoiner(System.lineSeparator());
        // print header
        outputJoiner.add("structure: '" + structure.getProteinIdentifier().getFullName() + "'")
                .add("chains: " + structure.chainsWithAminoAcids()
                        .map(Chain::getChainIdentifier)
                        .map(ChainIdentifier::getChainId)
                        .collect(Collectors.toList()))
                .add("total residues: " + structure.aminoAcids().count())
                .add("chain,res,aa,sse,energy,egor,sse_size,loop_fraction,rasa,plip_local_contacts," +
                        "plip_local_hbonds,plip_local_hydrophobic,plip_local_backbone,plip_long_range_contacts," +
                        "plip_long_range_hbonds,plip_long_range_hydrophobic,plip_long_range_backbone," +
                        "plip_betweenness,plip_closeness,plip_clusteringcoefficient,plip_hbonds_betweenness," +
                        "plip_hbonds_closeness,plip_hbonds_clusteringcoefficient,plip_hydrophobic_betweenness," +
                        "plip_hydrophobic_closeness,plip_hydrophobic_clusteringcoefficient,conv_betweenness," +
                        "conv_closeness,conv_clusteringcoefficient,plip_neighborhoods,conv_neighborhoods,prob,folds");

        structure.chainsWithAminoAcids()
                .forEach(chain -> {
                    List<String> output = structure.aminoAcids()
                            .map(aminoAcid -> {
                                boolean isProline = aminoAcid instanceof Proline;

                                SmoothedFeatureVector smoothedFeatureVector = aminoAcid.getFeature(SmoothedFeatureVector.class);
                                double loopFraction = aminoAcid.getFeature(LoopFraction.class).getLoopFraction();
                                Instance instance = createInstance(smoothedFeatureVector, loopFraction);
                                double prob = 0.0;
                                if (!isProline) {
                                    try {
                                        prob = model.distributionForInstance(normalize(instance))[0];
                                    } catch (Exception e) {
                                        e.printStackTrace();
                                    }
                                }

                                StringJoiner lineJoiner = new StringJoiner(",");
                                lineJoiner.add(aminoAcid.getParentChain().getChainIdentifier().getChainId())
                                        .add(aminoAcid.getResidueIdentifier().toString())
                                        .add(aminoAcid.getOneLetterCode())
                                        .add(aminoAcid.getFeature(GenericSecondaryStructure.class).getSecondaryStructure().getReducedRepresentation());
                                for (int i = 0; i < instance.numAttributes() - 1; i++) {
                                    lineJoiner.add(StandardFormat.format(instance.value(i)));
                                }
                                lineJoiner.add(StandardFormat.format(prob));
                                return lineJoiner.toString();
                            })
                            .sorted(Comparator.comparingDouble((String line) -> Double.valueOf(line.split(",")[line.split(",").length - 1])).reversed())
                            .collect(Collectors.toList());

                    int numberOfEarlyFoldingResidues = (int) (0.15 * (int) chain.aminoAcids().count());
                    int counter = 0;
                    for(int i = 0; i < chain.aminoAcids().count(); i++) {
                        outputJoiner.add(output.get(i) + "," + (counter < numberOfEarlyFoldingResidues ? "early" : "late"));
                        counter++;
                    }
                });

        // write output
        System.out.println("writing output to " + outputPath);
        Files.write(outputPath, outputJoiner.toString().getBytes());
    }

    public static void main(String[] args) throws IOException {
        if(args.length != 2) {
            System.out.println("incorrect number of arguments\n" +
                    "usage: java -jar efr.jar input output\n" +
                    "input: either '/a/path/to/some/file.pdb' or '1pdb' a pdb-id\n" +
                    "output: either '/a/absolute/path/output.csv' or 'relative.csv' a relative path");
            return;
        }

        String input = args[0];
        String output = args[1];
        Path outputPath = Paths.get(output);
        System.out.println();

        // parse structure - classic flow-control by exceptions
        Structure structure;
        try {
            ProteinIdentifier pdbId = IdentifierFactory.createProteinIdentifier(input);
            structure = StructureParser.fromProteinIdentifier(pdbId).parse();
            System.out.println("parsing structure by pdb-id: '" + pdbId + "'");
        } catch (Exception e) {
            System.out.println("parsing structure from file at: '" + input + "'");
            Path inputPath = Paths.get(input);
            structure = StructureParser.fromPath(inputPath).parse();
        }
        System.out.println();

        EarlyFoldingClassificationRunner instance = getInstance();
        instance.process(structure, outputPath);
    }

    private static final double[][] NORMALIZATION = new double[][] {
            { -10.847, 4.548 }, // energy
            { -12.017, 3.598 }, // egor
            { 8.991, 5.076 }, // sse_size
            { 0.43, 0.332 }, // loop_fraction
            { 0.308, 0.118 }, // rasa
            { 2.547, 1.282 }, // plip_local_contacts
            { 2.198, 1.197 }, // plip_local_hbonds
            { 0.2, 0.218 }, // plip_local_hydrophobic
            { 1.924, 1.191 }, // plip_local_backbone
            { 1.564, 1.05 }, // plip_long_range_contacts
            { 0.888, 0.878 }, // plip_long_range_hbonds
            { 0.56, 0.347 }, // plip_long_range_hydrophobic
            { 0.744, 0.851 }, // plip_long_range_backbone
            { 0.053, 0.024 }, // plip_betweenness
            { 0.21, 0.035 }, // plip_closeness
            { 0.185, 0.106 }, // plip_clusteringcoefficient
            { 0.075, 0.041 }, // plip_hbonds_betweenness
            { 0.142, 0.034 }, // plip_hbonds_closeness
            { 0.157, 0.116 }, // plip_hbonds_clusteringcoefficient
            { 0.072, 0.043 }, // plip_hydrophobic_betweenness
            { 0.151, 0.033 }, // plip_hydrophobic_closeness
            { 0.04, 0.065 }, // plip_hydrophobic_clusteringcoefficient
            { 0.041, 0.019 }, // conv_betweenness
            { 0.286, 0.054 }, // conv_closeness
            { 0.574, 0.059 }, // conv_clusteringcoefficient
            { 0.685, 0.309 }, // plip_neighborhoods
            { 1.927, 0.846 } // conv_neighborhoods
    };

    private Instance normalize(Instance instance) {
        for(int i = 0; i < instance.numAttributes() - 1; i++) {
            double[] normalization = NORMALIZATION[i];
            instance.setValue(i, (instance.value(i) - normalization[0]) / normalization[1]);
        }
        return instance;
    }

    private Instance createInstance(SmoothedFeatureVector smoothedFeatureVector, double loopFraction) {
        ArrayList<Attribute> attributes = new ArrayList<>();
        attributes.add(new Attribute("energy"));
        attributes.add(new Attribute("egor"));
        attributes.add(new Attribute("sse_size"));
        attributes.add(new Attribute("loop_fraction"));
        attributes.add(new Attribute("rasa"));
        attributes.add(new Attribute("plip_local_contacts"));
        attributes.add(new Attribute("plip_local_hbonds"));
        attributes.add(new Attribute("plip_local_hydrophobic"));
        attributes.add(new Attribute("plip_local_backbone"));
        attributes.add(new Attribute("plip_long_range_contacts"));
        attributes.add(new Attribute("plip_long_range_hbonds"));
        attributes.add(new Attribute("plip_long_range_hydrophobic"));
        attributes.add(new Attribute("plip_long_range_backbone"));
        attributes.add(new Attribute("plip_betweenness"));
        attributes.add(new Attribute("plip_closeness"));
        attributes.add(new Attribute("plip_clusteringcoefficient"));
        attributes.add(new Attribute("plip_hbonds_betweenness"));
        attributes.add(new Attribute("plip_hbonds_closeness"));
        attributes.add(new Attribute("plip_hbonds_clusteringcoefficient"));
        attributes.add(new Attribute("plip_hydrophobic_betweenness"));
        attributes.add(new Attribute("plip_hydrophobic_closeness"));
        attributes.add(new Attribute("plip_hydrophobic_clusteringcoefficient"));
        attributes.add(new Attribute("conv_betweenness"));
        attributes.add(new Attribute("conv_closeness"));
        attributes.add(new Attribute("conv_clusteringcoefficient"));
        attributes.add(new Attribute("plip_neighborhoods"));
        attributes.add(new Attribute("conv_neighborhoods"));
        attributes.add(new Attribute("folds"));
        Instances dataset = new Instances("classify", attributes, 1);
        dataset.setClassIndex(attributes.size() - 1);

        SparseInstance sparseInstance = new SparseInstance(28);
        sparseInstance.setValue(0, round(smoothedFeatureVector.getEnergy()));
        sparseInstance.setValue(1, round(smoothedFeatureVector.getEgor()));

        sparseInstance.setValue(2, round(smoothedFeatureVector.getSecondaryStructureElementSize()));
        sparseInstance.setValue(3, round(loopFraction));

        sparseInstance.setValue(4, round(smoothedFeatureVector.getRasa()));

        sparseInstance.setValue(5, round(smoothedFeatureVector.getLocalInteractions()));
        sparseInstance.setValue(6, round(smoothedFeatureVector.getLocalHydrogen()));
        sparseInstance.setValue(7, round(smoothedFeatureVector.getLocalHydrophobic()));
        sparseInstance.setValue(8, round(smoothedFeatureVector.getLocalBackbone()));

        sparseInstance.setValue(9, round(smoothedFeatureVector.getNonLocalInteractions()));
        sparseInstance.setValue(10, round(smoothedFeatureVector.getNonLocalHydrogen()));
        sparseInstance.setValue(11, round(smoothedFeatureVector.getNonLocalHydrophobic()));
        sparseInstance.setValue(12, round(smoothedFeatureVector.getNonLocalBackbone()));

        sparseInstance.setValue(13, round(smoothedFeatureVector.getBetweenness()));
        sparseInstance.setValue(14, round(smoothedFeatureVector.getCloseness()));
        sparseInstance.setValue(15, round(smoothedFeatureVector.getClusteringCoefficient()));

        sparseInstance.setValue(16, round(smoothedFeatureVector.getHydrogenBetweenness()));
        sparseInstance.setValue(17, round(smoothedFeatureVector.getHydrogenCloseness()));
        sparseInstance.setValue(18, round(smoothedFeatureVector.getHydrogenClusteringCoefficient()));

        sparseInstance.setValue(19, round(smoothedFeatureVector.getHydrophobicBetweenness()));
        sparseInstance.setValue(20, round(smoothedFeatureVector.getHydrophobicCloseness()));
        sparseInstance.setValue(21, round(smoothedFeatureVector.getHydrogenClusteringCoefficient()));

        sparseInstance.setValue(22, round(smoothedFeatureVector.getConvBetweenness()));
        sparseInstance.setValue(23, round(smoothedFeatureVector.getConvCloseness()));
        sparseInstance.setValue(24, round(smoothedFeatureVector.getConvClusteringCoefficient()));

        sparseInstance.setValue(25, round(smoothedFeatureVector.getDistinctNeighborhoods()));
        sparseInstance.setValue(26, round(smoothedFeatureVector.getConvDistinctNeighborhoods()));

        dataset.add(sparseInstance);
        sparseInstance.setDataset(dataset);
        sparseInstance.setClassMissing();

        return sparseInstance;
    }

    private double round(double value) {
        return (double) Math.round(value * 10000d) / 10000d;
    }
}
