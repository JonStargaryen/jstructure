package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.parser.ParsingException;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.nodes.Element;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.net.URL;
import java.util.List;
import java.util.Map;
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
    private static final GroupPrototypeParser INSTANCE = new GroupPrototypeParser();
    /**
     * Collection of already parsed group information.
     */
    private final Map<String, GroupPrototype> prototypes;

    private GroupPrototypeParser() {
        // load prototype files of all amino acids
        String basePath = "prototype/";
        this.prototypes = Stream.of("ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "HOH", "ILE", "LEU", "LYS", "MET",
                "MSE", "PHE", "PRO", "PYL", "SEC", "SER", "THR", "TRP", "TYR", "UNK", "VAL")
                .map(threeLetterCode -> basePath + threeLetterCode + ".xml")
                .map(this::getResourceAsStream)
                .map(this::createPrototype)
                .collect(Collectors.toConcurrentMap(GroupPrototype::getId,
                        Function.identity()));
    }

    private GroupPrototype createPrototype(InputStream inputStream) {
        try {
            try (InputStreamReader inputStreamReader = new InputStreamReader(inputStream)) {
                try (BufferedReader bufferedReader = new BufferedReader(inputStreamReader)) {
                    Document document = Jsoup.parse(bufferedReader.lines()
                            .collect(Collectors.joining(System.lineSeparator())));

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

                    String threeLetterCode = document.getElementsByTag("PDBx:three_letter_code")
                            .first()
                            .text();

                    List<Atom> prototypeAtoms = document.getElementsByTag("PDBx:chem_comp_atomCategory")
                            .first()
                            .children()
                            .stream()
                            .map(this::mapToPrototypeAtom)
                            .collect(Collectors.toList());

                    return new GroupPrototype(id,
                            name,
                            polymerType,
                            parentCompound,
                            oneLetterCode,
                            threeLetterCode,
                            prototypeAtoms);
                }
            }
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
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
        throw new ParsingException(polymerType + " is an unknown polymer type");
    }

    private GroupPrototype createPrototype(String id) {
        try {
            logger.info("creating prototype '{}' from online definition file", id);
            return createPrototype(new URL(String.format(DEFINITION_FETCH_URL, id)).openStream());
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    private InputStream getResourceAsStream(String filename) {
        return Thread.currentThread().getContextClassLoader().getResourceAsStream(filename);
    }

    public static GroupPrototypeParser getInstance() {
        return INSTANCE;
    }

    public GroupPrototype getPrototype(String id) {
        return prototypes.computeIfAbsent(id, this::createPrototype);
    }
}

