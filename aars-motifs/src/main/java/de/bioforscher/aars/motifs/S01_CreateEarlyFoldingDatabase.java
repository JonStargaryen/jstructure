package de.bioforscher.aars.motifs;

import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.testutil.FileUtils;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

public class S01_CreateEarlyFoldingDatabase {
    private static final Path BASE_DIR = Paths.get("/home/bittrich/git/phd_sb_repo/data/aars/EFoldMine_code/pdb-motifs/");

    public static void main(String[] args) {
        StringJoiner efoldmineScript = new StringJoiner(System.lineSeparator());
        List<String> lines = new BufferedReader(new InputStreamReader(Thread.currentThread().getContextClassLoader().getResourceAsStream("nrpdb.list"))).lines().collect(Collectors.toList());

        Structure lastStructure = null;
        String lastPdbId = null;

        for (String line : lines) {
            String[] split = line.split("\t");
            String pdbId = split[0];
            String chainId = split[1];
            String id = pdbId + "_" + chainId;
            System.out.println(id);

            try {
                Structure structure;
                if (pdbId.equals(lastPdbId)) {
                    structure = lastStructure;
                } else {
                    structure = StructureParser.fromPdbId(pdbId).parse();
                }

                Chain chain = structure.select().chainId(chainId).asChain();
                String sequence = chain.getAminoAcidSequence();

                FileUtils.write(BASE_DIR.resolve("sequences").resolve(id + ".fasta"),
                        ">" + id + System.lineSeparator() + sequence);
                efoldmineScript.add("python2 EFoldMine.py pdb-motifs/sequences/" + id + ".fasta -o pdb-motifs/ef-predictions/" + id + ".out");

                lastPdbId = pdbId;
                lastStructure = structure;
            } catch (Exception e) {
                System.err.println("could not process: " + line);
                e.printStackTrace();
            }
        }

        FileUtils.write(BASE_DIR.getParent().resolve("efoldmine-pdb-motifs.sh"), efoldmineScript.toString());
    }
}
