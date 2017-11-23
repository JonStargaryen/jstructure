package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.model.identifier.ChainIdentifier;
import de.bioforscher.jstructure.model.identifier.IdentifierFactory;
import de.bioforscher.jstructure.model.identifier.ProteinIdentifier;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.model.structure.nucleotide.Nucleotide;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.time.LocalDate;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.zip.GZIPInputStream;

/**
 * A minimalistic parser for structures in <tt>PDB</tt> format.
 * Created by S on 29.09.2016.
 */
public class StructureParser {
    private static final Logger logger = LoggerFactory.getLogger(StructureParser.class);
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
    private boolean minimalParsing;

    private Structure protein;
    private StringBuilder titleString;
    private Chain currentChain;
    private Group currentGroup;
    private List<Chain> terminatedChains;
    private boolean passedFirstModel;

    private StructureParser(OptionalSteps builder) {
        skipModels = builder.skipModels;
        strictMode = builder.strictMode;
        skipHydrogens = builder.skipHydrogenAtoms;
        minimalParsing = builder.minimalParsing;

        protein = new Structure(ProteinIdentifier.UNKNOWN_PROTEIN_IDENTIFIER);
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

                    // if id is missing, use hinted fall back
                    if(idIsMissing()) {
                        protein.setProteinIdentifier(builder.hintProteinName);
                    }
                    protein.setTitle(titleString.length() > 0 ? titleString.toString() : DEFAULT_PROTEIN_TITLE);
                }
            }

            if(builder.forceProteinName != null) {
                protein.setProteinIdentifier(builder.forceProteinName);
            }
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    private boolean idIsMissing() {
        return protein.getProteinIdentifier() == null || protein.getProteinIdentifier().getPdbId() == null ||
                protein.getProteinIdentifier().getPdbId().isEmpty() || protein.getProteinIdentifier().equals(ProteinIdentifier.UNKNOWN_PROTEIN_IDENTIFIER);
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

    private static final List<String> NUMBERS = IntStream.range(1, 10)
            .mapToObj(String::valueOf)
            .collect(Collectors.toList());

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
            // try to parse header line components - implicitly fallback values are provided by Structure's constructor
            try {
                String classification = line.substring(10, 50).trim();
                if(classification.isEmpty()) {
                    classification = ProteinIdentifier.UNKNOWN_PROTEIN_IDENTIFIER.getAdditionalName();
                }
                protein.setClassification(classification);
            } catch (Exception e) {
                logger.warn("failed to parse classification from line '{}'", line, e);
            }

            try {
                LocalDate depositionDate = LocalDate.parse(line.substring(50, 59).trim(),
                        StandardFormat.getPdbDateFormatInstance());
                if(depositionDate.isAfter(LocalDate.now())) {
                    depositionDate = depositionDate.minusYears(100);
                }
                protein.setDepositionDate(depositionDate);
            } catch (Exception e) {
                // potential legacy header: 'HEADER    MEMBRANE PROTEIN, TRANSPORT PROTEIN, SIG 20-JUN-16   5KK2             '
                try {
                    LocalDate depositionDate = LocalDate.parse(line.substring(51, 60).trim(),
                            StandardFormat.getPdbDateFormatInstance());
                    if(depositionDate.isAfter(LocalDate.now())) {
                        depositionDate = depositionDate.minusYears(100);
                    }
                    protein.setDepositionDate(depositionDate);
                } catch (Exception e2) {
                    logger.warn("failed to parse depositionDate from line '{}'", line, e2);
                }
            }

            try {
                ProteinIdentifier proteinIdentifier = IdentifierFactory.createProteinIdentifier(line.substring(62, 66));
                protein.setProteinIdentifier(proteinIdentifier);
            } catch (Exception e) {
                // potential legacy header: 'HEADER    MEMBRANE PROTEIN, TRANSPORT PROTEIN, SIG 20-JUN-16   5KK2
                try {
                    ProteinIdentifier proteinIdentifier = IdentifierFactory.createProteinIdentifier(line.substring(63, 67));
                    protein.setProteinIdentifier(proteinIdentifier);
                } catch (Exception e2) {
                    logger.warn("failed to parse the pdbId from line '{}'", line, e2);
                }
            }
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
            String atomName = line.substring(12, 16).trim();
            String pdbName = line.substring(17, 20).trim();
//            String elementString = "";
//            try {
//                elementString = line.substring(76, 78).trim();
//            } catch (StringIndexOutOfBoundsException e) {
//                // happens for malformed files like from CASP
//                logger.debug("missing element definition in line:{}{}",
//                        System.lineSeparator(),
//                        line);
//            }
//            if(elementString.isEmpty()) {
//                if(strictMode) {
//                    throw new ParsingException("PDB parsing failed for line:" + System.lineSeparator() + "'" + line + "'");
//                } else {
//                    logger.debug("PDB parsing failed for line:{}'{}'", System.lineSeparator(), line);
////                    elementString = Element.X.name();
            String elementString = atomName.substring(0,  1);

            // hydrogen atoms start with a digit: take 2nd position in that case
            if(NUMBERS.contains(elementString)) {
                elementString = atomName.substring(1, 2);
            }
//                }
//            }

            Element element = Element.valueOfIgnoreCase(elementString);
            if(skipHydrogens && element.isHydrogen()) {
                return;
            }

            String alternativeLocationIndicator = line.substring(16, 17).trim();
            String rawChainId = line.substring(21, 22);
            rawChainId = rawChainId.equals(" ") ? Chain.UNKNOWN_CHAIN.getChainIdentifier().getChainId() : rawChainId;
            ChainIdentifier chainId = IdentifierFactory.createChainIdentifier(protein.getProteinIdentifier(), rawChainId);
            int resNum = Integer.parseInt(line.substring(22, 26).trim());
            String insertionCode = line.substring(26, 27).trim();

            if(currentChain == null || !currentChain.getChainIdentifier().equals(chainId)) {
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

            if(currentGroup == null || currentGroup.getResidueIdentifier().getResidueNumber() != resNum ||
                    !currentGroup.getResidueIdentifier().getInsertionCode().equals(insertionCode) ||
                    !currentGroup.getParentChain().getChainIdentifier().equals(chainId)) {
                    // residue changed - create new group object and set reference
                    currentGroup = createGroup(pdbName, IdentifierFactory.createResidueIdentifier(resNum, insertionCode),
                            terminatedChains.contains(currentChain), minimalParsing);
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

            // we append the current group with additional atoms
            Atom atom = Atom.builder(element,
                    new double[] { Double.valueOf(line.substring(30, 38).trim()),
                            Double.valueOf(line.substring(38, 46).trim()),
                            Double.valueOf(line.substring(46, 54).trim())
                    })
                    .name(atomName)
                    .pdbSerial(Integer.valueOf(line.substring(6, 11).trim()))
                    .occupancy(occupancy)
                    .bfactor(bfactor)
                    .alternativeLocation(alternativeLocationIndicator)
                    .build();

            // 17/05/22 - stopping to skip alternative positions
            currentGroup.addAtom(atom);
        }

        if(line.startsWith(END_MODEL_PREFIX)) {
            //TODO handling of multiple models
            passedFirstModel = true;
            logger.debug("skipping models for {}", protein.getProteinIdentifier().getFullName());
        }
    }

    private static Group createGroup(String pdbName, ResidueIdentifier residueIdentifier, boolean ligand, boolean minimalParsing) {
        GroupPrototypeParser groupPrototypeParser = minimalParsing ? GroupPrototypeParser.getFastInstance() : GroupPrototypeParser.getInstance();
        GroupPrototype prototype = groupPrototypeParser.getPrototype(pdbName);

        // it is an amino acid
        if(prototype.getPolymerType() == GroupPrototype.PolymerType.PEPTIDE_LINKING || prototype.getPolymerType() ==
                GroupPrototype.PolymerType.PEPTIDE_LIKE) {
            return AminoAcid.Family.createAminoAcid(pdbName, residueIdentifier, ligand);
        }

        // it is a nucleotide
        //TODO potentially more polymer types for nucleotides
        if(prototype.getPolymerType() == GroupPrototype.PolymerType.NA_LINKING) {
            return Nucleotide.Family.createNucleotide(pdbName, residueIdentifier, ligand);
        }

        if(pdbName.equals(Water.THREE_LETTER_CODE)) {
            return new Water(residueIdentifier);
        }

        // it is neither
        Group group = new Group(prototype, residueIdentifier, ligand);
        // force the name parsed from the file, otherwise upon unknown ligands (and, thus, missing prototypes) this
        // information would get lost
        group.setThreeLetterCode(pdbName);
        return group;
    }

    Structure getProtein() {
        return protein;
    }

    public static OptionalSteps source(InputStream inputStream) {
        return new OptionalSteps(inputStream);
    }

    //TODO users could expect this to be a filepath too - change method name? decide internally on-the-fly?
    public static OptionalSteps source(String pdbId) {
        return new OptionalSteps(pdbId).hintProteinName(IdentifierFactory.createProteinIdentifier(pdbId));
    }

    public static OptionalSteps source(Path path) {
        return new OptionalSteps(path).hintProteinName(IdentifierFactory.createProteinIdentifier("", path.toFile().getName()));
    }

    public static class OptionalSteps {
        InputStream inputStream;
        private String pdbId;
        private Path path;
        boolean skipModels = true;
        boolean strictMode = false;
        boolean skipHydrogenAtoms = true;
        boolean approximateMissingAtoms = false;
        boolean minimalParsing = false;
        ProteinIdentifier forceProteinName;
        ProteinIdentifier hintProteinName;
        //TODO remove or place in global config
        private static Path localPdbDirectory = Paths.get("/var/local/pdb/");
        private static final int DEFAULT_CACHE_SIZE = 1000;

        public static void setLocalPdbDirectory(Path localPdbDirectory) {
            OptionalSteps.localPdbDirectory = localPdbDirectory;
        }

        public static Path getLocalPdbDirectory() {
            return localPdbDirectory;
        }

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

        public OptionalSteps minimalParsing(boolean minimalParsing) {
            this.minimalParsing = minimalParsing;
            return this;
        }

        public OptionalSteps strictMode(boolean strictMode) {
            this.strictMode = strictMode;
            return this;
        }

        public OptionalSteps skipHydrogenAtoms(boolean skipHydrogenAtoms) {
            this.skipHydrogenAtoms = skipHydrogenAtoms;
            return this;
        }

        public OptionalSteps forceProteinName(ProteinIdentifier proteinName) {
            this.forceProteinName = proteinName;
            return this;
        }

        public OptionalSteps hintProteinName(ProteinIdentifier proteinName) {
            this.hintProteinName = proteinName;
            return this;
        }

        public OptionalSteps approximateMissingAtoms(boolean approximateMissingAtoms) {
            this.approximateMissingAtoms = approximateMissingAtoms;
            //TODO impl
            return this;
        }

        public OptionalSteps cachedMode() {
            return cachedMode(DEFAULT_CACHE_SIZE);
        }

        public OptionalSteps cachedMode(int cacheSize) {
            //TODO impl: keep track of loaded entries, so 2nd invocation will return a map entry
            return this;
        }

        public Structure parse() {
            try {
                if (pdbId != null) {
                    Path pdbDirectory = OptionalSteps.localPdbDirectory;
                    try {
                        if (pdbDirectory != null && Files.isDirectory(pdbDirectory)) {
                            logger.debug("using local PDB to provide {}", pdbId);
                            // use local PDB if setup
                            pdbId = pdbId.toLowerCase();
                            String middle = pdbId.substring(1, 3);
                            this.inputStream = new GZIPInputStream(Files.newInputStream(pdbDirectory.resolve(middle).resolve("pdb" + pdbId + ".ent.gz")));
                        }
                    } catch (IOException e) {
//                        logger.warn("did not find file {} in local PDB at {}", pdbId, pdbDirectory, e);
                    }
                    if(this.inputStream == null) {
                        // no local PDB found - fetch from www
                        this.inputStream = new URL(String.format(PDB_FETCH_URL, pdbId)).openStream();
                    }
                }

                if (path != null) {
                    this.inputStream = Files.newInputStream(path);

                    if(forceProteinName == null) {
                        String fileName = path.toFile().getName();
                        // will cause file names containing multiple '.' to drop information: pdbFile.getName().split("\\.")[0]
                        int end = fileName.lastIndexOf(".");
                        end = end != -1 ? end : fileName.length();
                        forceProteinName = IdentifierFactory.createProteinIdentifier("", fileName.substring(0, end));
                    }
                }
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }

            return new StructureParser(this).getProtein();
        }
    }
}