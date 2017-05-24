package de.bioforscher.jstructure.model.structure.prototype;

import de.bioforscher.jstructure.model.structure.Atom;

import java.util.List;
import java.util.Optional;

/**
 * Created by bittrich on 5/24/17.
 */
public class GroupPrototype {
    private final String id;
    private final String name;
    private final Group.PolymerType polymerType;
    private final GroupPrototype parentCompound;
    private final String oneLetterCode;
    private final String threeLetterCode;
    private final List<Atom> prototypeAtoms;

    GroupPrototype(String id) {
        this.id = id;
    }

    public String getId() {
        return id;
    }

    public String getName() {
        return name;
    }

    public Group.PolymerType getPolymerType() {
        return polymerType;
    }

    public Optional<GroupPrototype> getParentCompound() {
        return Optional.ofNullable(parentCompound);
    }

    public Optional<String> getOneLetterCode() {
        return Optional.ofNullable(oneLetterCode);
    }

    public String getThreeLetterCode() {
        return threeLetterCode;
    }

    public List<Atom> getPrototypeAtoms() {
        return prototypeAtoms;
    }
}
