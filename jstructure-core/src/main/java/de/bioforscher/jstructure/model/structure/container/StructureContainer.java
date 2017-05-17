package de.bioforscher.jstructure.model.structure.container;

import de.bioforscher.jstructure.model.Identifiable;
import de.bioforscher.jstructure.model.feature.Featureable;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.structure.AtomRecordWriter;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;


/**
 * Marks this object as part of a hierarchy. Specifies no methods directly, but is useful so
 * {@link FeatureProvider} can be generically tied to classes implementing this
 * interface. As a rule of thumb, there are concrete groupings of elements (such as {@link Group} or {@link Chain} and
 * rather loose collections of similar elements such as {@link AtomContainer} which does not necessarily assumes any
 * clearAtoms connection between a set of atoms.
 * Created by S on 02.10.2016.
 */
public interface StructureContainer extends Featureable, AtomRecordWriter, Identifiable {

}
