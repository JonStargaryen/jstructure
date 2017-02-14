package de.bioforscher.aggregation;

import de.bioforscher.Constants;
import de.bioforscher.jstructure.model.structure.Protein;
import de.bioforscher.jstructure.model.structure.selection.Selection;
import de.bioforscher.jstructure.parser.ProteinParser;

import java.nio.file.Paths;

/**
 * Write selected non-redundant chains.
 * Not used as KinkFinder needs helix records.
 * Created by bittrich on 2/13/17.
 */
@Deprecated
public class S02_ExtractChains {
    public static void main(String[] args) {
        Constants.lines(Paths.get(Constants.PDBTM_ALPHA_NR_LIST))
                .filter(Constants.isCommentLine.negate())
                // some lines do not contain chain annotations
                .filter(line -> line.split("_").length == 2)
                // skip already processed lines
                .filter(line -> !Paths.get(Constants.CHAIN_PATH + line + Constants.PDB_SUFFIX).toFile().exists())
                // download files and write them to directory
                .forEach(line -> {
                    System.out.println(String.format("processing '%s'", line));
                    String pdbId = line.split("_")[0];
                    String chainId = line.split("_")[1];
                    Protein protein = ProteinParser.parsePDBFile(Paths.get(Constants.STRUCTURE_PATH + pdbId + Constants.PDB_SUFFIX));
                    byte[] output = Selection.on(protein)
                            .chainName(chainId)
                            .asChainContainer()
                            .composePDBRecord()
                            .getBytes();
                    Constants.write(Paths.get(Constants.CHAIN_PATH + line + Constants.PDB_SUFFIX), output);
                });
    }
}
