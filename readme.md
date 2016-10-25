# de.bioforscher.jstructure

*de.bioforscher.jstructure* is a light-weight module for structural bioinformatics. It provides
a stream-based API to work with protein structures in *PDB* format.

### Data Model
*de.bioforscher.jstructure* hierarchy is close to that BioJava of: a `Protein` contains any
 number of `Chain` objects which contain any number of `Residue` objects which
 contain any number of `Atom` objects. Each element of this hierarchy implements
 a `Container` interface, so it is possible to retrieve all registered children
 in a standardized manner. A `Residue` is an `AtomContainer`, so it provides
 a function `Stream<Atom> atoms()` which is used to get access to all `Atoms`
 linked to this particular `Residue`. A `Chain` is a `ResidueContainer` as well
 as an `AtomContainer`, thus one can retrieve all `Residue` objects by calling
 `Stream<Residue> residues()` and all `Atom` objects by `Stream<Atom> atoms()`.
  
 For more fine-grained selections of children elements, several convenience
 functions are implemented respectively specified by the `Container` interfaces.
 Each `AtomContainer` provides also access to methods such as 
 `Stream<Atom> backboneAtoms()` or `Stream<Atom> nonHydrogenAtoms()` which will
 employ custom `Predicate<Atom>` implementations on the raw atom stream provided
 by `Stream<Atom> atoms()`. `Stream<Pair<Atom>> atomPairsInContact(double distanceCutoff)`
 provides a convenient way to traverse all atom pairs in contact (defined by some
 distance threshold).
 
 Furthermore, each element is aware of its parent container. This reference is
 automatically set, when a children element (such as an `Atom`) is added to a
 parent container (in that case a `Residue`). This allows for individual atoms
 to return their *PDB* representation as an *ATOM* record (which depends
 information on the parent `Residue` as well as the parent `Chain`). In
 consequence every element of the hierarchy can compose its *PDB* representation
 by gathering all *ATOM* records of all `Atom` objects associated to it. This
 functionality is specified by the `AtomRecordWriter` interface.
 
 Each of these classes implements `FeatureContainer`, so arbitrary data can be
 attached to any one of them.

### Design Approach
* methods could primarily run on and return streams
* stream-returning methods are named by their content, 
 i.e. `Stream<Atom> atoms()` in contrast to `List<Atom> getAtoms()`
 
### Using The API
Coming soon!
