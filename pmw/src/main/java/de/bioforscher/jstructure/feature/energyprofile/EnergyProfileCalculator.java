package de.bioforscher.jstructure.feature.energyprofile;

import de.bioforscher.jstructure.mathematics.LinearAlgebra3D;
import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.container.GroupContainer;
import de.bioforscher.jstructure.model.structure.scheme.BetaCarbonRepresentationScheme;
import de.bioforscher.jstructure.model.structure.scheme.RepresentationScheme;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static de.bioforscher.jstructure.feature.energyprofile.EnergyProfileCalculator.FeatureNames.SOLVATION_ENERGY;

/**
 * The standard implementation of Frank Dressel's energy profiles for globular and trans-membrane proteins. A
 * characteristic energy value is assigned to each residue within the structure. This value describes the interactions
 * of this residue with its neighboring residues as well as the solvent.<br />
 * The computation is based on a unique preference of each amino acid to be either exposed to the solvent or be embedded
 * in the hydrophobic core of the protein. For each residue of the structure all neighboring residues are extracted
 * (i.e. they are less than 8 A apart) and their preferences are summed up.<br /><br />
 *
 * References:<br />
 * <ul>
 *     <li>1</li>
 *     <li>2</li>
 *     <li>3</li>
 * </ul>
 *
 * TODO implement routine for trans-membrane proteins
 * TODO implement interaction energy term
 * Created by bittrich on 12/15/16.
 */
public class EnergyProfileCalculator implements FeatureProvider {
    private static final Logger logger = LoggerFactory.getLogger(EnergyProfileCalculator.class);
    /**
     * 2 residues are considered to be in contact, when the euclidean distance of their beta-carbons is below 8.0 A.
     */
    private static final double DEFAULT_INTERACTION_CUTOFF = 8.0;
    private static final String BASE_PATH = "energyprofile/";
    private static final String GLOBULAR_SOLVATION_PATH = BASE_PATH + "globular-solvation.dat";
    private static Map<String, Double> globularSolvationData;

    public EnergyProfileCalculator() {
        if(globularSolvationData == null) {
            initializeLibrary();
        }
    }

    public enum FeatureNames {
        // the energy term describing interactions with the solvent
        SOLVATION_ENERGY
    }

    @Override
    public void process(Protein protein) {
//        processBySelectionAPI(protein); - ~4x slower than naive impl
        processNaively(protein);
    }

    void processNaively(Protein protein) {
        final RepresentationScheme representationScheme = new BetaCarbonRepresentationScheme();
        final double squaredDistanceCutoff = DEFAULT_INTERACTION_CUTOFF * DEFAULT_INTERACTION_CUTOFF;

        for(Group currentGroup : protein.getGroups()) {
            if(!currentGroup.isAminoAcid()) {
                continue;
            }

            double[] currentGroupCoordinates = representationScheme.determineRepresentingAtom(currentGroup).getCoordinates();
            double solvation = 0;
            double currentGroupSolvationValue = globularSolvationData.get(currentGroup.getThreeLetterCode());

            for(Group surroundingGroup : protein.getGroups()) {
                if(!surroundingGroup.isAminoAcid()) {
                    continue;
                }

                double[] surroundingGroupCoordinates = representationScheme.determineRepresentingAtom(surroundingGroup).getCoordinates();

                if(LinearAlgebra3D.distanceFast(currentGroupCoordinates, surroundingGroupCoordinates) > squaredDistanceCutoff) {
                    continue;
                }

                solvation += currentGroupSolvationValue + globularSolvationData.get(surroundingGroup.getThreeLetterCode());
            }

            currentGroup.setFeature(SOLVATION_ENERGY, solvation);
        }
    }

    void processBySelectionAPI(Protein protein) {
        // extract all amino acids
        GroupContainer residues = Selection.on(protein)
                .aminoAcids()
                .asGroupContainer();

        // extract all interacting residue pairs
        List<Pair<Group, Group>> interactingResiduePairs = Selection.pairsOn(residues)
                .betaCarbonDistance(DEFAULT_INTERACTION_CUTOFF)
                .asFilteredGroupPairs()
                .collect(Collectors.toList());

        // compute energy for each
        residues.groups().parallel().forEach(currentResidue -> {
            double solvationEnergy = interactingResiduePairs.stream()
                    // filter for entries describing interactions with this residue
                    .filter(pair -> pair.contains(currentResidue))
                    .flatMap(pair -> Stream.of(pair.getLeft(), pair.getRight()))
                    // retrieve amino acid instance for the given residue
                    .map(Group::getThreeLetterCode)
                    // safety net for mal-formed pdb names of amino acids
                    .map(String::toUpperCase)
                    // look up preference
                    .mapToDouble(globularSolvationData::get)
                    // sum overall neighboring residues
                    .sum();

            // assign value to currently processed residue
            currentResidue.setFeature(SOLVATION_ENERGY, solvationEnergy);
        });
    }

    private synchronized void initializeLibrary() {
        // parse globular solvation data
        globularSolvationData = new BufferedReader(new InputStreamReader(getResourceAsStream(GLOBULAR_SOLVATION_PATH))).lines()
                // skip header line
                .filter(line -> !line.startsWith("amino_acid"))
                .map(line -> line.split(" "))
                .collect(Collectors.toMap(key -> key[0], value -> -Math.log(Double.valueOf(value[1]) / Double.valueOf(value[2]))));
    }

    private InputStream getResourceAsStream(String filepath) {
        return Objects.requireNonNull(Thread.currentThread().getContextClassLoader().getResourceAsStream(filepath),
                "failed to findAny resource as InputStream");
    }
}
