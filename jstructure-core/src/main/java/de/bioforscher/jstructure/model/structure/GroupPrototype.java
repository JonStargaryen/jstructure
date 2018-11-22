package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;

import java.util.List;
import java.util.Optional;
import java.util.stream.Stream;

/**
 * Represents common properties of groups such as names, polymer-type, and a set of prototypical atom coordinates.
 * Created by bittrich on 5/24/17.
 */
public class GroupPrototype {
    public enum PolymerType {
        // non-polymer ligands
        NON_POLYMER,
        // amino acids
        PEPTIDE_LINKING,
        // nucleotides
        NA_LINKING,
        // the case for e.g. STA
        PEPTIDE_LIKE,
        PEPTIDE_TERMINUS,
        NA_TERMINUS,
        SACCHARIDE
    }

    /**
     * Gutteridge grouping (doi:10.1016/j.tibs.2005.09.006)
     */
    public enum GutteridgeGrouping {
        NONE,
        THIOL,
        HYDROXYL,
        GUANIDINIUM,
        IMIDAZOLE,
        AMINO,
        CARBOXYLATE,
        AMIDE
    }

    private final String id;
    private final String name;
    private final PolymerType polymerType;
    private final String parentCompound;
    private final String oneLetterCode;
    private final String threeLetterCode;
    private final List<Atom> prototypeAtoms;
    private final List<String> aromaticAtoms;
    private final List<GroupPrototypeParser.Bond> prototypeBonds;

    GroupPrototype(String id,
                   String name,
                   PolymerType polymerType,
                   String parentCompound,
                   String oneLetterCode,
                   String threeLetterCode,
                   List<Atom> prototypeAtoms,
                   List<String> aromaticAtomNames,
                   List<GroupPrototypeParser.Bond> prototypeBonds) {
        this.id = id;
        this.name = name;
        this.polymerType = polymerType;
        this.parentCompound = parentCompound;
        this.oneLetterCode = oneLetterCode;
        this.threeLetterCode = threeLetterCode;
        this.prototypeAtoms = prototypeAtoms;
        this.aromaticAtoms = aromaticAtomNames;
        this.prototypeBonds = prototypeBonds;
    }

    public String getId() {
        return id;
    }

    public String getName() {
        return name;
    }

    public PolymerType getPolymerType() {
        return polymerType;
    }

    public Optional<String> getParentCompound() {
        return Optional.ofNullable(parentCompound);
    }

    public Optional<String> getOneLetterCode() {
        // we resolve this here as it is a really rare case and would otherwise infer with the instantiation of the prototypes

        // the trivial case
        if(oneLetterCode != null) {
            return Optional.of(oneLetterCode);
        }

        // we may use the parent compound
        if(parentCompound != null) {
            GroupPrototype parentPrototype = GroupPrototypeParser.getInstance().getPrototype(parentCompound);
            return parentPrototype.getOneLetterCode();
        }

        return Optional.empty();
    }

    public String getThreeLetterCode() {
        return threeLetterCode;
    }

    public List<Atom> getPrototypeAtoms() {
        return prototypeAtoms;
    }

    public List<String> getAromaticAtoms() {
        return aromaticAtoms;
    }

    public List<GroupPrototypeParser.Bond> getPrototypeBonds() {
        return prototypeBonds;
    }

    public Stream<String> covalentNeighbors(String atomName) {
        return prototypeBonds.stream()
                .map(bond -> bond.getAtomName(atomName))
                .filter(Optional::isPresent)
                .map(Optional::get);
    }


    public double getMaximumAccessibleSurfaceArea() {
        return AminoAcid.Family.resolveGroupPrototype(this).getMaximumAccessibleSurfaceArea();
    }

    public GutteridgeGrouping getGutteridgeGrouping() {
        return AminoAcid.Family.resolveGroupPrototype(this).getGutteridgeGrouping();
    }

    public double getIsoelectricPoint() {
        return AminoAcid.Family.resolveGroupPrototype(this).getIsoelectricPoint();
    }

    @Override
    public String toString() {
        return "GroupPrototype{" +
                "id='" + id + '\'' +
                ", name='" + name + '\'' +
                ", polymerType=" + polymerType +
                ", parentCompound='" + parentCompound + '\'' +
                ", oneLetterCode='" + oneLetterCode + '\'' +
                ", threeLetterCode='" + threeLetterCode + '\'' +
                ", prototypeAtoms=" + prototypeAtoms +
                '}';
    }
}
