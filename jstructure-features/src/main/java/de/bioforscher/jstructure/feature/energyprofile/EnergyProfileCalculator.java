package de.bioforscher.jstructure.feature.energyprofile;

import de.bioforscher.jstructure.mathematics.LinearAlgebra;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;

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
 *     <li>F. Heinke and D. Labudde. MembraneContainer protein stability analyses by means of protein energy profiles in case of nephrogenic diabetes insipidus. Comput Math Methods Med, 2012:790281, February 2012.</li>
 *     <li>F. Heinke, A. Tuukkanen, and D. Labudde. Analysis of MembraneContainer Structure Stability in Diabetes insipidus. InTech, 2011, Edited by Kyuzi Kamoi,ISBN: 978-953-307-367-5,DOI: 10.5772/22258</li>
 *     <li>F. Heinke and D. Labudde. Predicting functionality of the non-expressed putative human OHCU decarboxylase by means of novel protein energy profile-based methods. Conference proceedings of 13. Nachwuchswissenschaftlerkonferenz, April 2012.</li>
 *     <li>F. Heinke, D. Stockmann, S. Schildbach, M. Langer and D. Labudde. eProS - A Bioinformatics Knowledgebase, Toolbox and Database for Characterizing Protein Function. BDAS 2015.</li>
 * </ul>
 *
 * TODO implement routine for trans-membrane proteins
 * Created by bittrich on 12/15/16.
 */
public class EnergyProfileCalculator extends FeatureProvider {
    private static final Logger logger = LoggerFactory.getLogger(EnergyProfileCalculator.class);
    /**
     * 2 residues are considered to be in contact, when the euclidean distance of their beta-carbons is below 8.0 A.
     */
    private static final double DEFAULT_INTERACTION_CUTOFF = 8.0;
    private static final String BASE_PATH = "energyprofile/";
    private static final String GLOBULAR_SOLVATION_PATH = BASE_PATH + "ep-solvation-globular.dat";
    private static Map<String, Double> globularSolvationData;

    public EnergyProfileCalculator() {
        if(globularSolvationData == null) {
            initializeLibrary();
        }
    }

    @Override
    protected void processInternally(Structure protein) {
        final double squaredDistanceCutoff = DEFAULT_INTERACTION_CUTOFF * DEFAULT_INTERACTION_CUTOFF;
        final List<AminoAcid> aminoAcids = protein.aminoAcids()
                .collect(Collectors.toList());

        for(AminoAcid currentGroup : aminoAcids) {
            Optional<Atom> currentGroupBetaCarbon = currentGroup.select()
                    .betaCarbonAtoms()
                    .asOptionalAtom();
            double[] currentGroupCoordinates = currentGroupBetaCarbon
                    .map(Atom::getCoordinates)
                    .orElseGet(() -> currentGroup.calculate().centroid().getValue());
            double solvation = 0;
            double currentGroupSolvationValue = resolve(globularSolvationData, currentGroup);

            for(AminoAcid surroundingGroup : aminoAcids) {
                if(currentGroup.equals(surroundingGroup)) {
                    continue;
                }

                Optional<Atom> surroundingGroupBetaCarbon = surroundingGroup.select()
                        .betaCarbonAtoms()
                        .asOptionalAtom();
                double[] surroundingGroupCoordinates = surroundingGroupBetaCarbon
                        .map(Atom::getCoordinates)
                        .orElseGet(() -> surroundingGroup.calculate().centroid().getValue());

                if(LinearAlgebra.on(currentGroupCoordinates).distanceFast(surroundingGroupCoordinates) > squaredDistanceCutoff) {
                    continue;
                }

                solvation += currentGroupSolvationValue + resolve(globularSolvationData, surroundingGroup);
            }

            currentGroup.getFeatureContainer().addFeature(new EnergyProfile(this, solvation));
        }
    }

    private double resolve(Map<String, Double> preferenceMap, AminoAcid group) {
        String threeLetterCode = group.getThreeLetterCode();
        // standard amino acid
        if(preferenceMap.containsKey(threeLetterCode)) {
            return preferenceMap.get(threeLetterCode);
        }

        // modified/non-standard - move to fallback - try annotated parent compound - else use plain ALA
        String fallback = group.getGroupPrototype()
                .getParentCompound()
                .orElse("ALA");
        logger.debug("encountered non-standard amino acid {}, using {} as fallback", threeLetterCode, fallback);

        return preferenceMap.getOrDefault(fallback, preferenceMap.get("ALA"));
    }

    private synchronized void initializeLibrary() {
        // parse globular solvation data
        globularSolvationData = getResourceAsStream(GLOBULAR_SOLVATION_PATH)
                // skip header line
                .filter(line -> !line.startsWith("amino_acid"))
                .map(line -> line.split(" "))
                .collect(Collectors.toMap(key -> key[0],
                        value -> -Math.log(Double.valueOf(value[1]) / Double.valueOf(value[2]))));
    }
}
