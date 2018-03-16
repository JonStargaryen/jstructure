package de.bioforscher.amyloid;

import de.bioforscher.jstructure.mathematics.IntegerInterval;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.StructureParser;

public class A01_ComparePrionStructures {
    private static final IntegerInterval RANGE = new IntegerInterval(218, 289);
    private static final String NATIVE = "2wvn";
    private static final String AMYLOID = "2kj3";

    public static void main(String[] args) {
        Chain nativeChain = StructureParser.fromPdbId(NATIVE)
                .parse()
                .getFirstChain()
                .select()
                .residueNumber(RANGE)
                .asIsolatedStructure()
                .getFirstChain();

        Chain amyloidChain = StructureParser.fromPdbId(AMYLOID)
                .parse()
                .getFirstChain()
                .select()
                .residueNumber(RANGE)
                .asIsolatedStructure()
                .getFirstChain();

        System.out.println(nativeChain);
        System.out.println(amyloidChain);
    }
}
