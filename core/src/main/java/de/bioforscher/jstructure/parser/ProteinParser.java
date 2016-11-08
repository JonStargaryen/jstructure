package de.bioforscher.jstructure.parser;

import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.Residue;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

/**
 * A minimalistic parser for protein structures in <tt>PDB</tt> format.
 * Created by S on 29.09.2016.
 */
public class ProteinParser {
    final Logger logger = LoggerFactory.getLogger(ProteinParser.class);
    private static final String HEADER_PREFIX = "HEADER";
    private static final String TITLE_PREFIX = "TITLE";
    private static final String ATOM_PREFIX = "ATOM";
    private static final String HETATM_PREFIX = "HETATM";
    private static final String DEFAULT_PROTEIN_TITLE = "NO DESCRIPTION";
    private static final String END_MODEL_PREFIX = "ENDMDL";
    /**
     * The PDB URL which can be used to fetch structures by ID (format this using the id and you are good to go).
     */
    private static final String PDB_FETCH_URL = "https://files.rcsb.org/download/%s.pdb";
    private Protein protein;
    private StringBuilder titleString;
    private Chain currentChain;
    private Residue currentResidue;
    private boolean passedFirstModel;

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
        return parsePDBFile(Files.newInputStream(pdbFile.toPath()), pdbFile.getName().split("\\.")[0]);
    }

    /**
     * Parses an {@link InputStream} of a <tt>PDB</tt> file. All other functions delegate to this method.
     * @param inputStream an input stream of <tt>PDB</tt> content
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
        currentResidue = null;

        // parse file
        new BufferedReader(new InputStreamReader(inputStream)).lines().forEach(this::parseLine);

        if(protein.getName() == null || protein.getName().isEmpty()) {
            protein.setName(proteinName);
        }
        protein.setTitle(titleString.length() > 0 ? titleString.toString() : DEFAULT_PROTEIN_TITLE);
        return protein;
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
                       .append(line.substring(10, 80).trim());
        }

        // parsing atom record - information we need is marked with an '*' - indirectly needed information (getChain/getResidue) marked with an '#'
        // some information will inform us about changing getChain/getResidue
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
            // we check for valid amino acids here, so no empty (nucleotide/ligand-only) chains are parsed
            String pdbName = line.substring(17, 20).trim();
//            System.out.println(line);
//            if(pdbName.length() < 3) {
//                // if pdbName is less than 3 chars long, this is no amino acid - e.g. 'U' vs 'THR'
//                // TODO at some point parsing nucleotides and ligands would be nice
//                return;
//            }

            // actually there is something to parse
            String chainId = line.substring(21, 22);
            int resNum = Integer.parseInt(line.substring(22, 26).trim());

            if(currentChain == null || !currentChain.getChainId().equals(chainId)) {
                // getChain changed - create new getChain object and set reference
                currentChain = new Chain(chainId);
                protein.addChain(currentChain);
            }

            if(currentResidue == null || currentResidue.getResidueNumber() != resNum) {
                // getResidue changed - create new getResidue object and set reference
                currentResidue = new Residue(pdbName, resNum);
                currentChain.addGroup(currentResidue);
            }

            // we append the current getResidue container with additional atoms
            Atom atom = new Atom(line.substring(12, 16).trim(),
                    Integer.valueOf(line.substring(6, 11).trim()),
                    line.substring(76, 78).trim(),
                    new double[] { Double.valueOf(line.substring(30, 38).trim()),
                            Double.valueOf(line.substring(38, 46).trim()),
                            Double.valueOf(line.substring(46, 54).trim())
                    },
                    Float.valueOf(line.substring(54, 60).trim()),
                    Float.valueOf(line.substring(60, 66).trim()));
            if(currentResidue.atoms().noneMatch(a -> a.getName().equals(atom.getName()))) {
                currentResidue.addAtom(atom);
            } else {
                logger.info("skipping atoms for {}", atom);
            }
        }

        if(line.startsWith(END_MODEL_PREFIX)) {
            passedFirstModel = true;
            logger.info("skipping models for {}", protein.getName());
        }
    }
}