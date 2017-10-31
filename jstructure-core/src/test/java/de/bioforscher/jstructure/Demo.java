package de.bioforscher.jstructure;

import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.aminoacid.Tryptophan;
import org.junit.Test;

import java.util.Arrays;
import java.util.Random;

public class Demo {
    @Test
    public void test() {
        // fetch a structure or load from local PDB if setup
        Structure structure = StructureParser.source("1brr").parse();

        // select a chain
        Chain chainB = structure.select()
                .chainId("B")
                .asChain();

        // or a residue
        AminoAcid aminoAcid1 = chainB.select()
                .residueNumber(60)
                .asAminoAcid();
        // and another one
        AminoAcid aminoAcid2 = chainB.select()
                .residueNumber(100)
                .asAminoAcid();

        // compute their distance
        System.out.println("distance of " + aminoAcid1 + " and " + aminoAcid2 + ": " +
                StandardFormat.format(aminoAcid1.calculate()
                        .centroid()
                        .distance(aminoAcid2.calculate()
                                .centroid())));

        // access amino acid-specific atoms
        chainB.select()
                .aminoAcids()
                .groupName("TRP")
                .asFilteredGroups()
                .map(Tryptophan.class::cast)
                .map(tryptophan -> tryptophan + " CG position: " +
                        Arrays.toString(tryptophan.getCg().getCoordinates()))
                .forEach(System.out::println);

        // compute features on-the-fly and resolve dependencies
        // e.g. assign some random value to each amino acid
        structure.aminoAcids()
                .forEach(aminoAcid -> aminoAcid.getFeatureContainer().addFeature(new Feature(new Random().nextDouble())));

        chainB.aminoAcids()
                .map(aminoAcid -> aminoAcid + " random feature: " +
                        StandardFormat.format(aminoAcid.getFeature(Feature.class).getValue()))
                .forEach(System.out::println);

        System.out.println("averages among chains:");
        structure.chainsWithAminoAcids()
                .map(chain -> chain.getChainIdentifier() + "'s average random feature: " +
                        StandardFormat.format(chain.aminoAcids()
                                .map(aminoAcid -> aminoAcid.getFeature(Feature.class))
                                .mapToDouble(Feature::getValue)
                                .average()
                                .getAsDouble()))
                .forEach(System.out::println);
    }

    static class Feature extends FeatureContainerEntry {
        private final double randomValue;

        public Feature(double randomValue) {
            // do not do this :p
            super(null);
            this.randomValue = randomValue;
        }

        public double getValue() {
            return randomValue;
        }
    }
}
