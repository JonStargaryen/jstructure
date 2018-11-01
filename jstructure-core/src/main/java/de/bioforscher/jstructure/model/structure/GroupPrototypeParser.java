package de.bioforscher.jstructure.model.structure;

import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.nodes.Element;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.net.URL;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Parses definition files locally stored for standard-cases or fetches them from the PDB.
 * Created by bittrich on 5/24/17.
 */
public class GroupPrototypeParser {
    private static final Logger logger = LoggerFactory.getLogger(GroupPrototypeParser.class);
    /**
     * The URL which can be used to fetch ligand information.
     */
    private static final String DEFINITION_FETCH_URL = "https://files.rcsb.org/ligands/view/%s.xml";
    /**
     * The regular instance.
     */
    private static final GroupPrototypeParser INSTANCE = new GroupPrototypeParser();
    /**
     * The 'fast' instance not struggling considering any ligands not known out-of-the-box.
     */
    private static final GroupPrototypeParser FAST_INSTANCE = new GroupPrototypeParser(true);
    /**
     * Collection of already parsed group information.
     */
    private final Map<String, GroupPrototype> prototypes;
    private boolean fastMode;
    private final GroupPrototype unknownLigand;

    /**
     * Used to determine which files are located at 'prototype/' as this cannot happen on-the-fly in bundled jars.
     * @param args
     */
    public static void main(String[] args) {

    }

    private GroupPrototypeParser() {
        this(false);
    }

    private GroupPrototypeParser(boolean fastMode) {
        String basePath = "prototype/";
        // load prototype files of all amino acids
        this.prototypes = Stream.of("ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU",
                "LYS", "MET", "MSE", "PHE", "PRO", "PYL", "SEC", "SER", "THR", "TRP", "TYR", "UNK", "VAL",
                // and nucleotides
                "A", "C", "DA", "DC", "DG", "G", "N", "T", "U",
                // and water and freaks
                "HOH", "UNL"
                )
                .map(threeLetterCode -> basePath + threeLetterCode + ".xml")
                .map(this::getResourceAsStream)
                .map(GroupPrototypeParser::getDocument)
                .map(this::createPrototype)
                .collect(Collectors.toConcurrentMap(GroupPrototype::getId,
                        Function.identity()));
        this.fastMode = fastMode;
        this.unknownLigand = prototypes.get("UNL");
    }

    private GroupPrototype createPrototype(Document document) {
        String id = document.getElementsByTag("PDBx:chem_comp")
                .first()
                .attr("id");

        String name = document.getElementsByTag("PDBx:name")
                .first()
                .text();

        GroupPrototype.PolymerType polymerType = document.getElementsByTag("PDBx:type")
                .stream()
                .map(Element::text)
                .map(text -> mapToPolymerType(text, id))
                .findFirst()
                .orElseThrow(() -> new ParsingException("definition file for '" + id + "' did not contain type tag"));

        String parentCompound = document.getElementsByTag("PDBx:mon_nstd_parent_comp_id")
                .stream()
                .map(Element::text)
                .findFirst()
                // parent compound may be not found for standard components
                .orElse(null);

        String oneLetterCode = document.getElementsByTag("PDBx:one_letter_code")
                .stream()
                .map(Element::text)
                .findFirst()
                // oneLetterCode will be undefined for 'real' ligands
                .orElse(null);

        String threeLetterCode = document.getElementsByTag("PDBx:chem_comp")
                .first()
                .attr("id");

        List<Atom> prototypeAtoms = new ArrayList<>();
        List<String> aromaticAtoms = new ArrayList<>();
        List<Bond> prototypeBonds = new ArrayList<>();

        try {
            prototypeAtoms = document.getElementsByTag("PDBx:chem_comp_atomCategory")
                    .first()
                    .children()
                    .stream()
                    .map(this::mapToPrototypeAtom)
                    .collect(Collectors.toList());

            aromaticAtoms = document.getElementsByTag("PDBx:chem_comp_atomCategory")
                    .first()
                    .children()
                    .stream()
                    .filter(element -> element.getElementsByTag("PDBx:pdbx_aromatic_flag").text().equals("Y"))
                    .map(element -> element.attr("atom_id"))
                    .collect(Collectors.toList());

            prototypeBonds = document.getElementsByTag("PDBx:chem_comp_bondCategory")
                    .first()
                    .children()
                    .stream()
                    .map(this::mapToPrototypeBond)
                    .collect(Collectors.toList());
        } catch (NullPointerException e) {
            // occurs for the unknown ligand ('UNL') which does not provide any prototype atoms
        }

        return new GroupPrototype(id,
                name,
                polymerType,
                parentCompound,
                oneLetterCode,
                threeLetterCode,
                prototypeAtoms,
                aromaticAtoms,
                prototypeBonds);
    }

    private Bond mapToPrototypeBond(Element element) {
        return new Bond(element.attr("atom_id_1"),
                element.attr("atom_id_2"),
                BondType.resolve(element.getElementsByTag("PDBx:value_order").text()));
    }

    public class Bond {
        private final String atomName1;
        private final String atomName2;
        private final BondType bondType;

        public Bond(String atomName1,
                    String atomName2,
                    BondType bondType) {
            this.atomName1 = atomName1;
            this.atomName2 = atomName2;
            this.bondType = bondType;
        }

        public String getAtomName1() {
            return atomName1;
        }

        public String getAtomName2() {
            return atomName2;
        }

        public BondType getBondType() {
            return bondType;
        }

        public boolean contains(String atomName) {
            return atomName1.equals(atomName) || atomName2.equals(atomName);
        }

        public Optional<String> getAtomName(String atomName) {
            if(atomName1.equals(atomName)) {
                return Optional.of(atomName2);
            } else if(atomName2.equals(atomName)) {
                return Optional.of(atomName1);
            } else {
                return Optional.empty();
            }
        }

        @Override
        public String toString() {
            return "Bond '" + atomName1 + "'-'" + atomName2 + "' : " + bondType;
        }
    }

    public enum BondType {
        SINGLE,
        DOUBLE,
        TRIPLE;

        public static BondType resolve(String description) {
            switch (description) {
                case "sing":
                    return BondType.SINGLE;
                case "doub":
                    return BondType.DOUBLE;
                case "trip":
                    return BondType.TRIPLE;
                default:
                    throw new UnsupportedOperationException("cannot handle bond type: " + description);
            }
        }
    }

    private Atom mapToPrototypeAtom(Element element) {
        return Atom.builder(de.bioforscher.jstructure.model.structure.Element.resolveElementSymbol(element.getElementsByTag("PDBx:type_symbol").text()),
                new double[]{
                        Double.valueOf(element.getElementsByTag("PDBx:pdbx_model_Cartn_x_ideal").text()),
                        Double.valueOf(element.getElementsByTag("PDBx:pdbx_model_Cartn_y_ideal").text()),
                        Double.valueOf(element.getElementsByTag("PDBx:pdbx_model_Cartn_z_ideal").text())
                })
                .name(element.getElementsByTag("PDBx:pdbx_component_atom_id").text())
                .pdbSerial(Integer.valueOf(element.getElementsByTag("PDBx:pdbx_ordinal").text()))
                .build();
    }

    private GroupPrototype.PolymerType mapToPolymerType(String polymerType, String id) {
        if(polymerType.contains("peptide linking") || (polymerType.contains("peptide") && polymerType.contains("linking"))) {
            return GroupPrototype.PolymerType.PEPTIDE_LINKING;
        }
        if(polymerType.contains("non-polymer")) {
            return GroupPrototype.PolymerType.NON_POLYMER;
        }
        if(polymerType.contains("NA linking")) {
            return GroupPrototype.PolymerType.NA_LINKING;
        }
        if(polymerType.contains("saccharide")) {
            return GroupPrototype.PolymerType.SACCHARIDE;
        }
        if(polymerType.contains("peptide-like")) {
            return GroupPrototype.PolymerType.PEPTIDE_LIKE;
        }
        if(polymerType.contains("peptide") && polymerType.contains("terminus")) {
            return GroupPrototype.PolymerType.PEPTIDE_TERMINUS;
        }
        throw new ParsingException("'" + polymerType + "' is an unknown polymer type for '" + id + "'");
    }

    private GroupPrototype createPrototype(String id) {
        return createPrototype(getDocument(id));
    }

    public static Document getDocument(String id) {
        try {
            logger.debug("creating prototype '{}' from online definition file", id);
            return getDocument(new URL(String.format(DEFINITION_FETCH_URL, id)).openStream());
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    private static Document getDocument(InputStream inputStream) {
        try {
            try (InputStreamReader inputStreamReader = new InputStreamReader(inputStream)) {
                try (BufferedReader bufferedReader = new BufferedReader(inputStreamReader)) {
                    return Jsoup.parse(bufferedReader.lines()
                            .collect(Collectors.joining(System.lineSeparator())));
                }
            }
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    private InputStream getResourceAsStream(String filepath) {
        return Objects.requireNonNull(Thread.currentThread().getContextClassLoader().getResourceAsStream(filepath),
                "failed to find resource as InputStream");
    }

    public static GroupPrototypeParser getInstance() {
        return INSTANCE;
    }

    public static GroupPrototypeParser getFastInstance() {
        return FAST_INSTANCE;
    }

    public GroupPrototype getPrototype(String id) {
        if(fastMode) {
            return prototypes.getOrDefault(id, unknownLigand);
        } else {
            // fetch and parse unknown ligand data if slow mode
            return prototypes.computeIfAbsent(id, this::createPrototype);
        }
    }
}

