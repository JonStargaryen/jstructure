package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.parser.ParsingException;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.nodes.Element;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Objects;
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
        // load prototype files of all amino acids
        String basePath = "prototype/";
        this.prototypes = Stream.of("ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "HOH", "ILE", "LEU",
                "LYS", "MET", "MSE", "PHE", "PRO", "PYL", "SEC", "SER", "THR", "TRP", "TYR", "UNK", "UNL", "VAL")
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
                .map(this::mapToPolymerType)
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

        List<Atom> prototypeAtoms;
        try {
            prototypeAtoms = document.getElementsByTag("PDBx:chem_comp_atomCategory")
                    .first()
                    .children()
                    .stream()
                    .map(this::mapToPrototypeAtom)
                    .collect(Collectors.toList());
        } catch (NullPointerException e) {
            // occurs for the unknown ligand ('UNL') which does not provide any prototype atoms
            prototypeAtoms = new ArrayList<>();
        }

        return new GroupPrototype(id,
                name,
                polymerType,
                parentCompound,
                oneLetterCode,
                threeLetterCode,
                prototypeAtoms);
    }

    private Atom mapToPrototypeAtom(Element element) {
        //TODO impl
        return null;
    }

    private GroupPrototype.PolymerType mapToPolymerType(String polymerType) {
        if(polymerType.contains("peptide linking")) {
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
        throw new ParsingException(polymerType + " is an unknown polymer type");
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

