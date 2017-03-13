package feature.topology;

import de.bioforscher.jstructure.feature.motif.SequenceMotif;
import de.bioforscher.jstructure.feature.motif.SequenceMotifAnnotator;
import de.bioforscher.jstructure.feature.topology.ANVIL;
import de.bioforscher.jstructure.feature.topology.Membrane;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import de.bioforscher.jstructure.parser.ProteinParser;
import org.junit.Before;
import org.junit.Test;

import java.util.List;

/**
 * Tests whether the potential calculation for membrane protein residues is working.
 * Created by S on 10.11.2016.
 */
public class MembranePotentialIntegrationTest {
    private Protein protein1brr;

    @Before
    public void setup() {
        protein1brr = ProteinParser.source("1brr").parse();
        new ANVIL().process(protein1brr);
        new SequenceMotifAnnotator().process(protein1brr);
    }

    @Test
    public void shouldComputeMembraneDistanceForEachResidue() {
        Membrane membrane = protein1brr.getFeature(Membrane.class, ANVIL.MEMBRANE);
        Selection.on(protein1brr)
                .aminoAcids()
                .asFilteredGroups()
                .map(residue -> residue.getParentChain().getChainId() + "-" +
                        residue.getThreeLetterCode() + "-" + residue.getResidueNumber() + " : " +
                        membrane.distanceToMembraneCenter(residue))
                .forEach(System.out::println);
    }

    @Test
    public void shouldComputeMembranePotentialForEachResidue() {
        Membrane membrane = protein1brr.getFeature(Membrane.class, ANVIL.MEMBRANE);
        Selection.on(protein1brr)
                .aminoAcids()
                .asFilteredGroups()
                .map(residue -> residue.getParentChain().getChainId() + "-" +
                        residue.getThreeLetterCode() + "-" + residue.getResidueNumber() + " : " +
                        membrane.computePotential(residue))
                .forEach(System.out::println);
    }

    @Test
    @SuppressWarnings("unchecked")
    public void shouldComputeMembranePotentialForEachSequenceMotif() {
        Membrane membrane = protein1brr.getFeature(Membrane.class, ANVIL.MEMBRANE);
        List<SequenceMotif> motifs = protein1brr.getFeature(List.class, SequenceMotifAnnotator.SEQUENCE_MOTIF);
        motifs.stream()
              .map(motif -> motif + " : " + membrane.computePotential(motif))
              .forEach(System.out::println);
    }
}
