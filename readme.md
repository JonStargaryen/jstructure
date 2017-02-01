# de.bioforscher.jstructure

*de.bioforscher.jstructure* is a light-weight module for structural bioinformatics. It provides
a stream-based API to work with protein structures in *PDB* format.

### Data Model
*de.bioforscher.jstructure* hierarchy is close to that BioJava of: a `Protein` contains any
 number of `Chain` objects which contain any number of `Group` objects which
 contain any number of `Atom` objects. Each element of this hierarchy implements
 a `Container` interface, so it is possible to retrieve all registered children
 in a standardized manner. A `Group` is an `AtomContainer`, so it provides
 a function `Stream<Atom> atoms()` which is used to get access to all `Atoms`
 linked to this particular `Group`. A `Chain` is a `GroupContainer` as well
 as an `AtomContainer`, thus one can retrieve all `Group` objects by calling
 `Stream<Group> groups()` and all `Atom` objects by `Stream<Atom> atoms()`.
  
 `Selection` provides a more fine-grained retrieval of child elements in a nature similar to
 that of a step-wise builder and will recognize the level of the selection.
 
 Furthermore, each element is aware of its parent container. This reference is
 automatically set, when a children element (such as an `Atom`) is added to a
 parent container (in that case a `Group`). This allows for individual atoms
 to return their *PDB* representation as an *ATOM* record (which depends
 information on the parent `Group` as well as the parent `Chain`). In
 consequence every element of the hierarchy can compose its *PDB* representation
 by gathering all *ATOM* records of all `Atom` objects associated to it. This
 functionality is specified by the `AtomRecordWriter` interface.
 
 Each of these classes implements `FeatureContainer`, so arbitrary data can be
 attached to any one of them.

### Feature Providers

These `AbstractFeatureProvider` implementations are tracked in the `FeatureProviderRegistry`.
This is useful respectively needed as the computation of some feature may require other
features to be present beforehand. E.g., in order to calculate the membrane topology, the 
accessible surface area has to be computed first. With this simple implementation of a 
service registry, every implementation can get an instance of the appropriate classes to
compute missing features on-the-fly without the the user having to remember all requirements
while also not directly tying `AbstractFeatureProvider` instances together, but rather connecting
them by the features they require and provide.

### Using The API
        // fetch/parse structure by id
        Protein protein = ProteinParser.parseProteinById("1brr");
        
        // print coordinates of all alanines
        Selection.on(protein)
                .aminoAcids(AminoAcidFamily.ALANINE)
                .asFilteredGroups()
                .map(Group::composePDBRecord)
                .forEach(System.out::println);
        
        // count all alanines
        double alanineRatio = Selection.on(protein)
                .aminoAcids(AminoAcidFamily.ALANINE)
                .asFilteredGroups()
                .count() / (double) protein.getSize() * 100.0;
        
        // store count
        protein.setFeature("ALANINE_RATIO", alanineRatio);
        
        // retrieve it
        System.out.printf("alanine ratio: %3.2f%%", protein.getFeatureAsDouble("ALANINE_RATIO"));
        
        
        // compute the ASA by a suitable provider
        FeatureProviderRegistry.resolve("ACCESSIBLE_SURFACE_AREA").process(protein);
        
        // print values
        protein.aminoAcids()
                .map(group -> group.getFeatureAsDouble("ACCESSIBLE_SURFACE_AREA"))
                .forEach(System.out::println);

### The design Module

The `design` package contains example classes on how to aggregate information and 
store respectively analyze it.

These classes are used to tweak the `core` package and spot shortcomings of the
current design. Also they provide examples on how to employ the capabilities of it.
Be aware, that everything in all modules is suspect to drastic changes as the library
grows and requirements change.