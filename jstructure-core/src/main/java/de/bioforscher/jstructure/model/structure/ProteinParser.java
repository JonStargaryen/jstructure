package de.bioforscher.jstructure.model.structure;

import de.bioforscher.jstructure.model.structure.aminoacid.*;
import de.bioforscher.jstructure.model.structure.identifier.ChainIdentifier;
import de.bioforscher.jstructure.model.structure.identifier.IdentifierFactory;
import de.bioforscher.jstructure.model.structure.identifier.ProteinIdentifier;
import de.bioforscher.jstructure.model.structure.identifier.ResidueIdentifier;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.zip.GZIPInputStream;

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
    private boolean minimalParsing;

    private Protein protein;
    private StringBuilder titleString;
    private Chain currentChain;
    private Group currentGroup;
    private List<Chain> terminatedChains;
    private boolean passedFirstModel;

    private ProteinParser(OptionalSteps builder) {
        skipModels = builder.skipModels;
        strictMode = builder.strictMode;
        skipHydrogens = builder.skipHydrogenAtoms;
        minimalParsing = builder.minimalParsing;

        protein = new Protein(ProteinIdentifier.UNKNOWN_PROTEIN_ID);
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

    private boolean idIsMissing() {
        return protein.getPdbId() == null || protein.getPdbId().getPdbId() == null ||
                protein.getPdbId().getPdbId().isEmpty() || protein.getPdbId().equals(ProteinIdentifier.UNKNOWN_PROTEIN_ID);
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
            protein.setPdbId(IdentifierFactory.createProteinIdentifier(line.substring(62, 66)));
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

            String alternativeLocationIndicator = line.substring(16, 17).trim();
            String pdbName = line.substring(17, 20).trim();
            ChainIdentifier chainId = IdentifierFactory.createChainIdentifier(protein.getPdbId(), line.substring(21, 22));
            int resNum = Integer.parseInt(line.substring(22, 26).trim());
            String insertionCode = line.substring(26, 27).trim();

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

            if(currentGroup == null || currentGroup.getResidueIdentifier().getResidueNumber() != resNum ||
                    !currentGroup.getResidueIdentifier().getInsertionCode().equals(insertionCode) ||
                    !currentGroup.getParentChain().getChainId().equals(chainId)) {
                    // residue changed - create new group object and set reference
                    currentGroup = createGroup(pdbName, IdentifierFactory.createResidueIdentifier(resNum, insertionCode),
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

            // we append the current group with additional atoms
            Atom atom = Atom.builder(element,
                    new double[] { Double.valueOf(line.substring(30, 38).trim()),
                            Double.valueOf(line.substring(38, 46).trim()),
                            Double.valueOf(line.substring(46, 54).trim())
                    })
                    .name(line.substring(12, 16).trim())
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
            logger.debug("skipping models for {}", protein.getPdbId().getFullName());
        }
    }

    private Group createGroup(String pdbName, ResidueIdentifier residueIdentifier, boolean ligand) {
        GroupPrototypeParser groupPrototypeParser = minimalParsing ? GroupPrototypeParser.getFastInstance() : GroupPrototypeParser.getInstance();
        GroupPrototype prototype = groupPrototypeParser.getPrototype(pdbName);

        // it is an amino acid
        if(prototype.getPolymerType() == GroupPrototype.PolymerType.PEPTIDE_LINKING || prototype.getPolymerType() ==
                GroupPrototype.PolymerType.PEPTIDE_LIKE)
        switch(pdbName.toUpperCase()) {
            case Alanine.THREE_LETTER_CODE:
                return new Alanine(residueIdentifier, ligand);
            case Arginine.THREE_LETTER_CODE:
                return new Arginine(residueIdentifier, ligand);
            case Asparagine.THREE_LETTER_CODE:
                return new Asparagine(residueIdentifier, ligand);
            case AsparticAcid.THREE_LETTER_CODE:
                return new AsparticAcid(residueIdentifier, ligand);
            case Cysteine.THREE_LETTER_CODE:
                return new Cysteine(residueIdentifier, ligand);
            case GlutamicAcid.THREE_LETTER_CODE:
                return new GlutamicAcid(residueIdentifier, ligand);
            case Glutamine.THREE_LETTER_CODE:
                return new Glutamine(residueIdentifier, ligand);
            case Glycine.THREE_LETTER_CODE:
                return new Glycine(residueIdentifier, ligand);
            case Histidine.THREE_LETTER_CODE:
                return new Histidine(residueIdentifier, ligand);
            case Isoleucine.THREE_LETTER_CODE:
                return new Isoleucine(residueIdentifier, ligand);
            case Leucine.THREE_LETTER_CODE:
                return new Leucine(residueIdentifier, ligand);
            case Lysine.THREE_LETTER_CODE:
                return new Lysine(residueIdentifier, ligand);
            case Methionine.THREE_LETTER_CODE:
                return new Methionine(residueIdentifier, ligand);
            case Phenylalanine.THREE_LETTER_CODE:
                return new Phenylalanine(residueIdentifier, ligand);
            case Proline.THREE_LETTER_CODE:
                return new Proline(residueIdentifier, ligand);
            case Pyrrolysine.THREE_LETTER_CODE:
                return new Pyrrolysine(residueIdentifier, ligand);
            case Selenocysteine.THREE_LETTER_CODE:
                return new Selenocysteine(residueIdentifier, ligand);
            case Selenomethionine.THREE_LETTER_CODE:
                return new Selenomethionine(residueIdentifier, ligand);
            case Serine.THREE_LETTER_CODE:
                return new Serine(residueIdentifier, ligand);
            case Threonine.THREE_LETTER_CODE:
                return new Threonine(residueIdentifier, ligand);
            case Tryptophan.THREE_LETTER_CODE:
                return new Tryptophan(residueIdentifier, ligand);
            case Tyrosine.THREE_LETTER_CODE:
                return new Tyrosine(residueIdentifier, ligand);
            case Valine.THREE_LETTER_CODE:
                return new Valine(residueIdentifier, ligand);
            default:
                return new UnknownAminoAcid(pdbName, residueIdentifier, ligand);
        }

        // it is a nucleotide
        if(prototype.getPolymerType() == GroupPrototype.PolymerType.NA_LINKING) {
            //TODO impl
        }

        // it is neither
        Group group = new Group(prototype, residueIdentifier, ligand);
        // force the name parsed from the file, otherwise upon unknown ligands (and, thus, missing prototypes) this
        // information would get lost
        group.setThreeLetterCode(pdbName);
        return group;
    }

    Protein getProtein() {
        return protein;
    }

    //TODO the local pdb approach handling is horrible
    //TODO traverse-all function for local pdb
    public static OptionalSteps localPdb(String pdbId) {
        String middle = pdbId.substring(1, 3);
        try {
            Path pdbDirectory = OptionalSteps.localPdbDirectory;
            if(pdbDirectory == null || !Files.isDirectory(pdbDirectory)) {
                throw new IllegalArgumentException("local pdb directory '" + pdbDirectory + "' is not set or not valid - use ProteinParser.OptionalSteps.setLocalPdbDirectory to set up");
            }
            InputStream inputStream = Files.newInputStream(pdbDirectory.resolve(middle).resolve("pdb" + pdbId + ".ent.gz"));
            return new OptionalSteps(new GZIPInputStream(inputStream)).hintProteinName(IdentifierFactory.createProteinIdentifier(pdbId));
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
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
        private static Path localPdbDirectory;
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
                        forceProteinName = IdentifierFactory.createProteinIdentifier("", fileName.substring(0, fileName.lastIndexOf(".")));
                    }
                }
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }

            return new ProteinParser(this).getProtein();
        }
    }
}