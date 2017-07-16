package de.bioforscher.jstructure.mutation.impl;

import de.bioforscher.jstructure.model.identifier.ChainIdentifier;
import de.bioforscher.jstructure.model.identifier.ResidueIdentifier;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import de.bioforscher.jstructure.mutation.MutatorService;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.stream.Collectors;

/**
 * Use scwrl4 to place the new side-chain reasonably.
 * @see <a href="http://dunbrack.fccc.edu/scwrl4/index.php">http://dunbrack.fccc.edu/scwrl4/index.php</a>
 * Created by bittrich on 7/14/17.
 */
public class ScwrlMutatorServiceImpl implements MutatorService {
    private static final Path SCWRL_DIRECTORY = Paths.get("/usr/local/bin/scwrl4/");
    private static final String SCWRL_COMMAND = SCWRL_DIRECTORY.resolve("Scwrl4").toFile().getAbsolutePath();

    @Override
    public Structure mutateAminoAcid(Structure originalProtein,
                                     ChainIdentifier chainIdentifier,
                                     ResidueIdentifier residueIdentifier,
                                     AminoAcid.Family targetAminoAcid) {
        try {
            AminoAcid aminoAcidToMutate = originalProtein.select()
                    .chainName(chainIdentifier.getChainId())
                    .residueNumber(residueIdentifier.getResidueNumber())
                    .asAminoAcid();

            Path tmpDirectory = Files.createTempDirectory("scwrl");

            // write structure file to tmp
            Path structurePath = tmpDirectory.resolve("original.pdb");
            Files.write(structurePath, originalProtein.getPdbRepresentation().getBytes());

            // write sequence file to tmp
            Path sequencePath = tmpDirectory.resolve("sequence.txt");
            Files.write(sequencePath, composeMutateScwrlSequence(originalProtein, aminoAcidToMutate, targetAminoAcid).getBytes());

            // expected output file
            Path outputPath = tmpDirectory.resolve("output.pdb");

            ProcessBuilder processBuilder = new ProcessBuilder(SCWRL_COMMAND,
                    "-i",
                    structurePath.toFile().getAbsolutePath(),
                    "-o",
                    outputPath.toFile().getAbsolutePath(),
                    "-s",
                    sequencePath.toFile().getAbsolutePath());

            processBuilder.start().waitFor();

            return StructureParser.source(outputPath).parse();
        } catch (Exception e) {
            //TODO error-handling
            throw new RuntimeException(e);
        }
    }

    /**
     * Composes a sequence in scwrl format, used to mutate a single side-chain.
     * <pre>
     -s <sequencefilename> [optional]
     This flag is followed by a sequence file. The sequence should have the same number of residues in it as the input
     backbone. White space, carriage returns, and numbers are ignored. Lower-case letters in the sequence indicate that
     the Cartesian coordinates for the corresponding residues are to be left untouched, and will be treated as steric
     boundaries only for the other side chains.

     Examples:
     SDERYCNM - full SCWRL side-chain replacement
     SdERYCNM - input residue (aspartate) is left where as is.
     SxERYCNM - input residue (aspartate) is left where as is.
     * </pre>
     *
     * @param originalProtein the original container
     * @param aminoAcidToMutate the original amino acid, providing the position to mutate
     * @param targetAminoAcid the target of this mutation
     * @return the sequence used as input for scwrl
     */
    public String composeMutateScwrlSequence(Structure originalProtein, AminoAcid aminoAcidToMutate, AminoAcid.Family targetAminoAcid) {
        return originalProtein.aminoAcids()
                .map(aminoAcid -> {
                    if(aminoAcidToMutate == aminoAcid) {
                        // uppercase amino acids in the sequence will be manipulated
                        return targetAminoAcid.getOneLetterCode();
                    } else {
                        // lowercase olc will be ignored
                        return aminoAcid.getOneLetterCode().toLowerCase();
                    }
                })
                .collect(Collectors.joining());
    }
}
