package de.bioforscher.jstructure.parser;

import de.bioforscher.jstructure.model.structure.*;
import de.bioforscher.jstructure.model.structure.identifier.PdbChainId;
import de.bioforscher.jstructure.model.structure.identifier.PdbId;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

/**
 * A minimalistic parser for protein structures in <tt>PDB</tt> format.
 * Created by S on 29.09.2016.
 */
public class ProteinParser {
    private static final Logger logger = LoggerFactory.getLogger(ProteinParser.class);
    private static final String HEADER_PREFIX = "HEADER";
    private static final String TITLE_PREFIX = "TITLE";
    private static final String DEFAULT_PROTEIN_TITLE = "NO DESCRIPTION";
    private static final String END_MODEL_PREFIX = "ENDMDL";
    private static final String TER_PREFIX = "TER";
    /**
     * The PDB URL which can be used to fetch structures by ID (format this using the id and you are good to go).
     */
    private static final String PDB_FETCH_URL = "https://files.rcsb.org/download/%s.pdb";

    private boolean skipModels;
    private boolean strictMode;
    private boolean skipHydrogens;

    private Protein protein;
    private StringBuilder titleString;
    private Chain currentChain;
    private Group currentGroup;
    private List<Chain> terminatedChains;
    private boolean passedFirstModel;

    private ProteinParser(OptionalSteps builder) {
        skipModels = builder.skipModels;
        strictMode = builder.strictMode;
        skipHydrogens = builder.skipHydrogens;

        protein = new Protein();
        // 'initialize' title field as it tends to be split over multiple lines - thus, we have to append previous results when we find further entries
        titleString = new StringBuilder();

        // null fields, so the same instance can be used multiple times
        currentChain = null;
        currentGroup = null;

        // keep track of processed TER records
        terminatedChains = new ArrayList<>();

        // parse file
        try {
            try(InputStreamReader inputStreamReader = new InputStreamReader(builder.inputStream)) {
                try(BufferedReader bufferedReader = new BufferedReader(inputStreamReader)) {
                    bufferedReader.lines().forEach(this::parseLineChecked);

                    if(protein.getPdbId() == null || protein.getPdbId().getPdbId() == null || protein.getPdbId().getPdbId().isEmpty()) {
                        protein.setPdbId(builder.hintProteinName);
                    }
                    protein.setTitle(titleString.length() > 0 ? titleString.toString() : DEFAULT_PROTEIN_TITLE);
                }
            }

            if(builder.forceProteinName != null) {
                protein.setPdbId(builder.forceProteinName);
            }
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    private void parseLineChecked(String line) {
        try {
            parseLine(line);
        } catch (StringIndexOutOfBoundsException | NumberFormatException e) {
            if(strictMode) {
                throw new ParsingException("PDB parsing failed for line:" + System.lineSeparator() + "'" + line + "'", e);
            } else {
                logger.debug("PDB parsing failed for line:{}'{}'{}cause: {}", System.lineSeparator(), line, System.lineSeparator(), e);
            }
        }
    }

    /**
     * Parses a single line of a <tt>PDB</tt> file.
     * @param line the line to processUniProtId
     */
    private void parseLine(String line) {
        //TODO this is kinda hacky, however that way only the first model is parsed in any case
        if(passedFirstModel) {
            return;
        }
        // indices taken from: ftp://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_Letter.pdf
        // their column definition has certain offset to the definition of String#substring(int, int)

        // found header record
        // 63 - 66       IDcode        idCode            This identifier is unique within the PDB.
        if(line.startsWith(HEADER_PREFIX)) {
            protein.setPdbId(PdbId.createFromPdbId(line.substring(62, 66)));
        }

        // handling title records
        // 11 - 80       String        title         Title of the experiment.
        if(line.startsWith(TITLE_PREFIX)) {
            // trim to omit tailing white-spaces
            // extra whitespace to ensure that words are separated
            // maybe some StringJoiner is the way to go
            titleString.append(titleString.length() == 0 ? "" : " ")
                    .append(line.substring(10, line.length() < 80 ? line.length() : 80).trim());
        }

        if(line.startsWith(TER_PREFIX)) {
            // mark chain as terminated - everything parsed from now on, associated to this chain will be an HETATM
            Chain chainToTerminate = protein.select()
                    .chainName(line.length() > 22 ? line.substring(21, 22) : "?")
                    .asOptionalChain()
                    .orElse(currentChain);
            terminatedChains.add(chainToTerminate);
        }

        // parsing atom record - information we need is marked with an '*' - indirectly needed information (chain/residue) marked with an '#'
        // some information will inform us about changing chain/residue
		/*	COLUMNS        DATA TYPE     FIELD        DEFINITION
			-------------------------------------------------------------------------------------
			1 - 6          Record name   "ATOM  "
		*	7 - 11   	   Integer       serial       Atom serial number.
		*	13 - 16        Atom          name         Atom name.
			17             Character     altLoc       Alternate location indicator.
		#	18 - 20        Residue name  resName      Residue name.
		#	22             Character     chainID      Chain identifier.
		#	23 - 26        Integer       resSeq       Residue sequence number.
			27             AChar         iCode        Code for insertion of residues.
		*	31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
		*	39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
		*	47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
			55 - 60        Real(6.2)    occupancy     Occupancy.
			61 - 66        Real(6.2)    tempFactor    Temperature factor.
		*	77 - 78        LString(2)   element       Element symbol, right justified.
			79 - 80        LString(2)   charge        Charge on the atom */
        if(line.startsWith(Atom.ATOM_PREFIX) || line.startsWith(Atom.HETATM_PREFIX)) {
            String elementString = line.substring(76, 78).trim();
            if(elementString.isEmpty()) {
                if(strictMode) {
                    throw new ParsingException("PDB parsing failed for line:" + System.lineSeparator() + "'" + line + "'");
                } else {
                    logger.debug("PDB parsing failed for line:{}'{}'", System.lineSeparator(), line);
                    elementString = Element.X.name();
                }
            }

            Element element = Element.valueOfIgnoreCase(elementString);
            if(skipHydrogens && element.isHydrogen()) {
                return;
            }

            String pdbName = line.substring(17, 20).trim();
            PdbChainId chainId = PdbChainId.createFromChainId(protein.getPdbId(), line.substring(21, 22));
            int resNum = Integer.parseInt(line.substring(22, 26).trim());

            if(currentChain == null || !currentChain.getChainId().equals(chainId)) {
                Optional<Chain> selectedChain = protein.select()
                        .chainName(chainId.getChainId())
                        .asOptionalChain();
                if(selectedChain.isPresent()) {
                    // chain already present - just an het-group not directly connected
                    currentChain = selectedChain.get();
                } else {
                    // chain changed - create new chain object and set reference
                    currentChain = new Chain(chainId);
                    protein.addChain(currentChain);
                }
            }

            if(currentGroup == null || currentGroup.getResidueNumber() != resNum ||
                    !currentGroup.getParentChain().getChainId().equals(chainId)) {
                // residue changed - create new group object and set reference
                currentGroup = new Group(pdbName,
                        resNum,
                        CIFParser.parseLigandInformation(pdbName),
                        terminatedChains.contains(currentChain));
                currentChain.addGroup(currentGroup);
            }

            float occupancy;
            try {
                occupancy = Float.valueOf(line.substring(54, 60).trim());
            } catch (NumberFormatException e) {
                if(strictMode) {
                    throw new ParsingException(e);
                } else {
                    logger.debug("missing occupancy in line{}'{}'", System.lineSeparator(), line);
                    occupancy = Atom.DEFAULT_OCCUPANCY;
                }
            }

            float bfactor;
            try {
                bfactor = Float.valueOf(line.substring(60, 66).trim());
            } catch (NumberFormatException e) {
                if(strictMode) {
                    throw new ParsingException(e);
                } else {
                    logger.debug("missing bfactor in line{}'{}'", System.lineSeparator(), line);
                    bfactor = Atom.DEFAULT_BFACTOR;
                }
            }

            // we append the current getResidue container with additional atoms
            Atom atom = new Atom(line.substring(12, 16).trim(),
                    Integer.valueOf(line.substring(6, 11).trim()),
                    element,
                    new double[] { Double.valueOf(line.substring(30, 38).trim()),
                            Double.valueOf(line.substring(38, 46).trim()),
                            Double.valueOf(line.substring(46, 54).trim())
                    },
                    occupancy,
                    bfactor);
            if(currentGroup.atoms().noneMatch(a -> a.getName().equals(atom.getName()))) {
                currentGroup.addAtom(atom);
            } else {
                //TODO find solution for this - maybe only the atom with the highest occupancy should be kept
                logger.debug("skipping atoms for {}", atom);
            }
        }

        if(line.startsWith(END_MODEL_PREFIX)) {
            //TODO handling of multiple models
            passedFirstModel = true;
            logger.debug("skipping models for {}", protein.getPdbId().getFullName());
        }
    }

    Protein getProtein() {
        return protein;
    }

    public static OptionalSteps source(InputStream inputStream) {
        return new OptionalSteps(inputStream);
    }

    //TODO users could expect this to be a filepath too - change method name? decide internally on-the-fly?
    public static OptionalSteps source(String pdbId) {
        return new OptionalSteps(pdbId).hintProteinName(PdbId.createFromPdbId(pdbId));
    }

    public static OptionalSteps source(Path path) {
        return new OptionalSteps(path).hintProteinName(PdbId.createFromName(path.toFile().getName()));
    }

    public static class OptionalSteps {
        InputStream inputStream;
        private String pdbId;
        private Path path;
        boolean skipModels = true;
        boolean strictMode = false;
        boolean skipHydrogens = true;
        PdbId forceProteinName;
        PdbId hintProteinName;
        private static final int DEFAULT_CACHE_SIZE = 1000;

        OptionalSteps(InputStream inputStream) {
            this.inputStream = inputStream;
        }

        OptionalSteps(String pdbId) {
            this.pdbId = pdbId;
        }

        OptionalSteps(Path path) {
            this.path = path;
        }

        public OptionalSteps skipModels(boolean skipModels) {
            this.skipModels = skipModels;
            return this;
        }

        public OptionalSteps strictMode(boolean strictMode) {
            this.strictMode = strictMode;
            return this;
        }

        public OptionalSteps skipHydrogens(boolean skipHydrogens) {
            this.skipHydrogens = skipHydrogens;
            return this;
        }

        public OptionalSteps forceProteinName(PdbId proteinName) {
            this.forceProteinName = proteinName;
            return this;
        }

        public OptionalSteps hintProteinName(PdbId proteinName) {
            this.hintProteinName = proteinName;
            return this;
        }

        public OptionalSteps cachedMode() {
            return cachedMode(DEFAULT_CACHE_SIZE);
        }

        public OptionalSteps cachedMode(int cacheSize) {
            //TODO impl: keep track of loaded entries, so 2nd invocation will return a map entry
            return this;
        }

        public Protein parse() {
            try {
                if (pdbId != null) {
                    this.inputStream = new URL(String.format(PDB_FETCH_URL, pdbId)).openStream();
                }

                if (path != null) {
                    this.inputStream = Files.newInputStream(path);

                    if(forceProteinName == null) {
                        String fileName = path.toFile().getName();
                        // will cause file names containing multiple '.' to drop information: pdbFile.getName().split("\\.")[0]
                        forceProteinName = PdbId.createFromName(fileName.substring(0, fileName.lastIndexOf(".")));
                    }
                }
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }

            return new ProteinParser(this).getProtein();
        }
    }
}