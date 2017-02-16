package de.bioforscher.jstructure.parser;

import de.bioforscher.jstructure.model.structure.*;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

/**
 * A minimalistic parser for protein structures in <tt>PDB</tt> format.
 * Created by S on 29.09.2016.
 */
public class ProteinParser {
    private final Logger logger = LoggerFactory.getLogger(ProteinParser.class);
    private static final String HEADER_PREFIX = "HEADER";
    private static final String TITLE_PREFIX = "TITLE";
    public static final String ATOM_PREFIX = "ATOM  ";
    public static final String HETATM_PREFIX = "HETATM";
    private static final String DEFAULT_PROTEIN_TITLE = "NO DESCRIPTION";
    private static final String END_MODEL_PREFIX = "ENDMDL";
    private static final String TER_PREFIX = "TER";
    /**
     * The PDB URL which can be used to fetch structures by ID (format this using the id and you are good to go).
     */
    private static final String PDB_FETCH_URL = "https://files.rcsb.org/download/%s.pdb";
    private Protein protein;
    private StringBuilder titleString;
    private Chain currentChain;
    private Group currentGroup;
    private boolean passedFirstModel;
    private List<Chain> terminatedChains;
    //TODO move to more mature options
    private static boolean skipHydrogens = true;

    /**
     * Fetches a <tt>PDB</tt> file from the <tt>PDB</tt> based on a <tt>PDB</tt> id.
     * @param pdbId the 4 digit pdb code to fetch
     * @return the parsed protein
     * @see ProteinParser#PDB_FETCH_URL
     * @see ProteinParser#parsePDBFile(InputStream, String)
     */
    public static Protein parseProteinById(String pdbId) {
        try {
            return new ProteinParser().parsePDBFile(new URL(String.format(PDB_FETCH_URL, pdbId)).openStream(), pdbId);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    /**
     * Parses a local file based on a filepath.
     * @param path path of the file to parse
     * @return the parsed protein
     * @see ProteinParser#parsePDBFile(InputStream, String)
     */
    public static Protein parsePDBFile(Path path) {
        try {
            return new ProteinParser().parsePDBFile(path.toFile());
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    /**
     * Parses a local file based on a filepath.
     * @param filepath filepath of the file to parse
     * @return the parsed protein
     * @see ProteinParser#parsePDBFile(InputStream, String)
     */
    public static Protein parsePDBFile(String filepath) {
        return parsePDBFile(Paths.get(filepath));
    }

    /**
     * Parses a local file based on a file object.
     * @param pdbFile reference of the file to parse
     * @return the parsed protein
     * @see ProteinParser#parsePDBFile(InputStream, String)
     * @throws IOException when IO fails at some point
     */
    private Protein parsePDBFile(File pdbFile) throws IOException {
        try(InputStream inputStream = Files.newInputStream(pdbFile.toPath())) {
            String fileName = pdbFile.getName();
            return parsePDBFile(inputStream,
                    // will cause file names containing multiple '.' to drop information: pdbFile.getName().split("\\.")[0]
                    fileName.substring(0, fileName.lastIndexOf("."))
            );
        }
    }

    /**
     * Parses an {@link InputStream} of a <tt>PDB</tt> file. All other functions delegate to this method.
     * @param inputStream an input select of <tt>PDB</tt> content
     * @param proteinName the desired name - will only be assigned, when the parsed lines do not contain a valid
     *                    <tt>HEADER</tt> record
     * @return a {@link Protein} instance wrapping all (relevant) parsed information
     * @throws IOException when IO fails at some point
     */
    private Protein parsePDBFile(InputStream inputStream, String proteinName) throws IOException {
        protein = new Protein();
        // 'initialize' title field as it tends to be split over multiple lines - thus, we have to append previous results when we find further entries
        titleString = new StringBuilder();

        // null fields, so the same instance can be used multiple times
        currentChain = null;
        currentGroup = null;

        // keep track of processed TER records
        terminatedChains = new ArrayList<>();

        // parse file
        try(InputStreamReader inputStreamReader = new InputStreamReader(inputStream)) {
            try(BufferedReader bufferedReader = new BufferedReader(inputStreamReader)) {
                bufferedReader.lines().forEach(this::parseLineChecked);

                if(protein.getName() == null || protein.getName().isEmpty()) {
                    protein.setName(proteinName);
                }
                protein.setTitle(titleString.length() > 0 ? titleString.toString() : DEFAULT_PROTEIN_TITLE);
                return protein;
            }
        }
    }

    private void parseLineChecked(String line) {
        try {
            parseLine(line);
        } catch (StringIndexOutOfBoundsException | NumberFormatException e) {
            throw new ParsingException("PDB parsing failed for line:" + System.lineSeparator() + line, e);
        }
    }

    /**
     * Parses a single line of a <tt>PDB</tt> file.
     * @param line the line to process
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
            protein.setName(line.substring(62, 66));
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
            Chain chainToTerminate = Selection.on(protein)
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
        if(line.startsWith(ATOM_PREFIX) || line.startsWith(HETATM_PREFIX)) {
            String elementString = line.substring(76, 78).trim();
            if(elementString.isEmpty()) {
                throw new ParsingException("PDB parsing failed for line:" + System.lineSeparator() + line );
            }

            Element element = Element.valueOfIgnoreCase(elementString);
            if(skipHydrogens && element.isHydrogen()) {
                return;
            }

            String pdbName = line.substring(17, 20).trim();
            String chainId = line.substring(21, 22);
            int resNum = Integer.parseInt(line.substring(22, 26).trim());

            if(currentChain == null || !currentChain.getChainId().equals(chainId)) {
                Optional<Chain> selectedChain = Selection.on(protein)
                        .chainName(chainId)
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

            if(currentGroup == null || currentGroup.getResidueNumber() != resNum || !currentGroup.getParentChain().getChainId().equals(chainId)) {
                // residue changed - create new group object and set reference
                currentGroup = new Group(pdbName,
                        resNum,
                        CIFParser.parseLigandInformation(pdbName),
                        terminatedChains.contains(currentChain));
                currentChain.addGroup(currentGroup);
            }

            // we append the current getResidue container with additional atoms
            Atom atom = new Atom(line.substring(12, 16).trim(),
                    Integer.valueOf(line.substring(6, 11).trim()),
                    element,
                    new double[] { Double.valueOf(line.substring(30, 38).trim()),
                            Double.valueOf(line.substring(38, 46).trim()),
                            Double.valueOf(line.substring(46, 54).trim())
                    },
                    Float.valueOf(line.substring(54, 60).trim()),
                    Float.valueOf(line.substring(60, 66).trim()));
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
            logger.debug("skipping models for {}", protein.getName());
        }
    }
}