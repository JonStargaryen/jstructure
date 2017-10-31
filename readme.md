# de.bioforscher.jstructure

*de.bioforscher.jstructure* is a light-weight library for structural bioinformatics. It provides
a stream-based API to work with macromolecular structures in *PDB* format.

### Data Model
*de.bioforscher.jstructure* hierarchy is close to that BioJava of: a `Structure` contains any
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
 by gathering all *ATOM* records of all `Atom` objects associated to it.
 
 To any of these classes arbitrary data can be attached to.

### Using The API
        // fetch a structure or load from local PDB if setup
        Structure structure = StructureParser.source("1brr").parse();
        
        // select a chain
        Chain chainB = structure.select()
                .chainId("B")
                .asChain();

        // or a residue
        AminoAcid aminoAcid1 = chainB.select()
                .residueNumber(60)
                .asAminoAcid();
        // and another one
        AminoAcid aminoAcid2 = chainB.select()
                .residueNumber(100)
                .asAminoAcid();

        // compute their distance
        System.out.println("distance of " + aminoAcid1 + " and " + aminoAcid2 + ": " +
                StandardFormat.format(aminoAcid1.calculate()
                        .centroid()
                        .distance(aminoAcid2.calculate()
                                .centroid())));

        // access amino acid-specific atoms
        chainB.select()
                .aminoAcids()
                .groupName("TRP")
                .asFilteredGroups()
                .map(Tryptophan.class::cast)
                .map(tryptophan -> tryptophan + " CG position: " +
                        Arrays.toString(tryptophan.getCg().getCoordinates()))
                .forEach(System.out::println);

        // compute features on-the-fly and resolve dependencies
        // e.g. assign some random value to each amino acid
        structure.aminoAcids()
                .forEach(aminoAcid -> aminoAcid.getFeatureContainer().addFeature(new Feature(new Random().nextDouble())));

        chainB.aminoAcids()
                .map(aminoAcid -> aminoAcid + " random feature: " +
                        StandardFormat.format(aminoAcid.getFeature(Feature.class).getValue()))
                .forEach(System.out::println);

        System.out.println("averages among chains:");
        structure.chainsWithAminoAcids()
                .map(chain -> chain.getChainIdentifier() + "'s average random feature: " +
                        StandardFormat.format(chain.aminoAcids()
                                .map(aminoAcid -> aminoAcid.getFeature(Feature.class))
                                .mapToDouble(Feature::getValue)
                                .average()
                                .getAsDouble()))
                .forEach(System.out::println);

### Feature Providers

Several `FeatureProvider` implementations are provided which allow the computation or 
annotation of values such as secondary structure information, accessible surface area values,
 membrane topology, evolutionary information or UniProt data.
 
Dependencies between them are resolved automatically and the user can request features which
will be computed on-the-fly, when they are not already present.

### Alignments

...

### Design guidelines
* no public method takes `null` as argument, none will return `null`