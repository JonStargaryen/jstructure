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
* methods are `Stream`-first and avoid returning `Collections`
* stream-returning methods are named by their content, 
 i.e. `Stream<Atom> atoms()` in contrast to `List<Atom> getAtoms()`
 
### Using The API
    // fetch/parse structure by id
     Protein protein = ProteinParser.parseProteinById("1brr");
    
     // print all residue pairs whose distance is less than 8.0 A
     protein.residuePairsInContact(8.0)
            .forEach(System.out::println);
    
     // custom filters
     Predicate<Residue> alanineFilter = new AminoAcidFilter(Collections.singletonList(ALANINE));
     // print coordinates of all alanines
     protein.residues()
            .filter(alanineFilter)
            .map(Residue::composePDBRecord)
            .forEach(System.out::println);
    
     // count all alanines
     double alanineRatio = protein.residues()
             .filter(alanineFilter)
             .count() / (double) protein.getSize() * 100.0;
    
     // store count
     protein.setFeature(FeatureNames.ALANINE_RATIO, alanineRatio);
    
     // retrieve it
     System.out.printf("alanine ratio: %3.2f%%", protein.getDoubleFeature(FeatureNames.ALANINE_RATIO));

Examples of extending the core API are given in the `pmw` module which contains 
algorithms to predict or annotate certain features of a `Protein` such as the 
accessible surface area, sequence motifs, secondary structure elements, the 
of membrane proteins. Basically, any information that can be computed based on 
a protein sequence or structure. Furthermore, structure reconstruction algorithms
are provided which can generate an all-atom representation of a `Protein` merely 
by providing distance and constraints between particular amino acid positions in
a sequence.

These classes are used to tweak the `core` package and spot shortcomings of the
current design. Also they provide examples on how to employ the capabilities of it.

The `design` package contains example classes on how to aggregate information and 
store respectively analyze it.