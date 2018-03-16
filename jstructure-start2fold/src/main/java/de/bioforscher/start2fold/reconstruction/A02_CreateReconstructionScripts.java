package de.bioforscher.start2fold.reconstruction;

import de.bioforscher.start2fold.Start2FoldConstants;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

public class A02_CreateReconstructionScripts {
    private static final Path SCRIPT_DIRECTORY = Paths.get("/home/bittrich/git/phd_sb_repo/data/start2fold/reconstruction/");

    public static void main(String[] args) throws IOException {
        List<Path> maps = Files.list(Paths.get("/home/bittrich/git/phd_sb_repo/data/start2fold/reconstruction/maps/"))
                .collect(Collectors.toList());

        int bins = 8;
        for (int binIndex = 0; binIndex < bins; binIndex++) {
            StringJoiner stringJoiner = new StringJoiner(System.lineSeparator());
            for (int fileIndex = 0; fileIndex < maps.size(); fileIndex++) {
                if(fileIndex % bins != binIndex) {
                    continue;
                }

                Path mapPath = maps.get(fileIndex);
                String[] split = mapPath.toFile().getName().split("-");
                String stfId = split[0];

                String strategy = split[1];
                if(!strategy.equals("efr") && !strategy.equals("random")) {
                    continue;
                }
                // create directory
                stringJoiner.add("mkdir /tmp/reconstruction/");
                stringJoiner.add("mkdir /tmp/reconstruction/" + binIndex + "/");

                // create confold call
                stringJoiner.add("/home/bittrich/programs/confold_v1.0/confold.pl " +
                        "-rrtype ca " +
                        "-seq /home/bittrich/git/phd_sb_repo/data/start2fold/reconstruction/fasta/" + stfId + ".fasta " +
                        "-rr " + mapPath.toFile().getAbsolutePath() + " " +
                        "-ss /home/bittrich/git/phd_sb_repo/data/start2fold/reconstruction/sse/" + stfId + ".sse " +
                        "-o /tmp/reconstruction/" + binIndex + "/");

                String percentage = split[2];
                String samplingId = split[3].split("\\.")[0];
                // align to native structure
                // stfId, strategy, percentage, samplingId, confoldRunId
                stringJoiner.add("echo \"" + stfId + "," + strategy + "," + percentage + "," + samplingId + ",1\" >> confold-" + binIndex + ".out");
                stringJoiner.add("tmalign /tmp/reconstruction/" + binIndex + "/stage1/" + stfId + "_model1.pdb /home/bittrich/git/phd_sb_repo/data/start2fold/pdb/" + stfId + ".pdb | grep RMSD >> confold-" + binIndex + ".out");
                stringJoiner.add("echo \"" + stfId + "," + strategy + "," + percentage + "," + samplingId + ",2\" >> confold-" + binIndex + ".out");
                stringJoiner.add("tmalign /tmp/reconstruction/" + binIndex + "/stage1/" + stfId + "_model2.pdb /home/bittrich/git/phd_sb_repo/data/start2fold/pdb/" + stfId + ".pdb | grep RMSD >> confold-" + binIndex + ".out");
                stringJoiner.add("echo \"" + stfId + "," + strategy + "," + percentage + "," + samplingId + ",3\" >> confold-" + binIndex + ".out");
                stringJoiner.add("tmalign /tmp/reconstruction/" + binIndex + "/stage1/" + stfId + "_model3.pdb /home/bittrich/git/phd_sb_repo/data/start2fold/pdb/" + stfId + ".pdb | grep RMSD >> confold-" + binIndex + ".out");
                stringJoiner.add("echo \"" + stfId + "," + strategy + "," + percentage + "," + samplingId + ",4\" >> confold-" + binIndex + ".out");
                stringJoiner.add("tmalign /tmp/reconstruction/" + binIndex + "/stage1/" + stfId + "_model4.pdb /home/bittrich/git/phd_sb_repo/data/start2fold/pdb/" + stfId + ".pdb | grep RMSD >> confold-" + binIndex + ".out");
                stringJoiner.add("echo \"" + stfId + "," + strategy + "," + percentage + "," + samplingId + ",5\" >> confold-" + binIndex + ".out");
                stringJoiner.add("tmalign /tmp/reconstruction/" + binIndex + "/stage1/" + stfId + "_model5.pdb /home/bittrich/git/phd_sb_repo/data/start2fold/pdb/" + stfId + ".pdb | grep RMSD >> confold-" + binIndex + ".out");

                // cleanup
                stringJoiner.add("rm -r /tmp/reconstruction/" + binIndex + "/");
            }
            Start2FoldConstants.write(SCRIPT_DIRECTORY.resolve("confold-" + binIndex + ".sh"),
                    stringJoiner.toString());
        }

        for (int binIndex = 0; binIndex < bins; binIndex++) {
            StringJoiner stringJoiner = new StringJoiner(System.lineSeparator());
            for (int fileIndex = 0; fileIndex < maps.size(); fileIndex++) {
                if(fileIndex % bins != binIndex) {
                    continue;
                }

                Path mapPath = maps.get(fileIndex);
                String[] split = mapPath.toFile().getName().split("-");
                String stfId = split[0];

                String strategy = split[1];
                if(!strategy.equals("residues")) {
                    continue;
                }

                // create directory
                stringJoiner.add("mkdir /tmp/reconstruction/");
                stringJoiner.add("mkdir /tmp/reconstruction/" + binIndex + "/");

                // create confold call
                stringJoiner.add("/home/bittrich/programs/confold_v1.0/confold.pl " +
                        "-rrtype ca " +
                        "-seq /home/bittrich/git/phd_sb_repo/data/start2fold/reconstruction/fasta/" + stfId + ".fasta " +
                        "-rr " + mapPath.toFile().getAbsolutePath() + " " +
                        "-ss /home/bittrich/git/phd_sb_repo/data/start2fold/reconstruction/sse/" + stfId + ".sse " +
                        "-o /tmp/reconstruction/" + binIndex + "/");

                String percentage = split[2];
                String samplingId = split[3].split("\\.")[0];
                // align to native structure
                // stfId, strategy, percentage, samplingId, confoldRunId
                stringJoiner.add("echo \"" + stfId + "," + strategy + "," + percentage + "," + samplingId + ",1\" >> confold-" + binIndex + "-residues.out");
                stringJoiner.add("tmalign /tmp/reconstruction/" + binIndex + "/stage1/" + stfId + "_model1.pdb /home/bittrich/git/phd_sb_repo/data/start2fold/pdb/" + stfId + ".pdb | grep RMSD >> confold-" + binIndex + "-residues.out");
                stringJoiner.add("echo \"" + stfId + "," + strategy + "," + percentage + "," + samplingId + ",2\" >> confold-" + binIndex + "-residues.out");
                stringJoiner.add("tmalign /tmp/reconstruction/" + binIndex + "/stage1/" + stfId + "_model2.pdb /home/bittrich/git/phd_sb_repo/data/start2fold/pdb/" + stfId + ".pdb | grep RMSD >> confold-" + binIndex + "-residues.out");
                stringJoiner.add("echo \"" + stfId + "," + strategy + "," + percentage + "," + samplingId + ",3\" >> confold-" + binIndex + "-residues.out");
                stringJoiner.add("tmalign /tmp/reconstruction/" + binIndex + "/stage1/" + stfId + "_model3.pdb /home/bittrich/git/phd_sb_repo/data/start2fold/pdb/" + stfId + ".pdb | grep RMSD >> confold-" + binIndex + "-residues.out");
                stringJoiner.add("echo \"" + stfId + "," + strategy + "," + percentage + "," + samplingId + ",4\" >> confold-" + binIndex + "-residues.out");
                stringJoiner.add("tmalign /tmp/reconstruction/" + binIndex + "/stage1/" + stfId + "_model4.pdb /home/bittrich/git/phd_sb_repo/data/start2fold/pdb/" + stfId + ".pdb | grep RMSD >> confold-" + binIndex + "-residues.out");
                stringJoiner.add("echo \"" + stfId + "," + strategy + "," + percentage + "," + samplingId + ",5\" >> confold-" + binIndex + "-residues.out");
                stringJoiner.add("tmalign /tmp/reconstruction/" + binIndex + "/stage1/" + stfId + "_model5.pdb /home/bittrich/git/phd_sb_repo/data/start2fold/pdb/" + stfId + ".pdb | grep RMSD >> confold-" + binIndex + "-residues.out");

                // cleanup
                stringJoiner.add("rm -r /tmp/reconstruction/" + binIndex + "/");
            }
            Start2FoldConstants.write(SCRIPT_DIRECTORY.resolve("confold-" + binIndex + "-residues.sh"),
                    stringJoiner.toString());
        }

        for (int binIndex = 0; binIndex < 2; binIndex++) {
            StringJoiner stringJoiner = new StringJoiner(System.lineSeparator());
            for (int fileIndex = 0; fileIndex < maps.size(); fileIndex++) {
                if(fileIndex % 2 != binIndex) {
                    continue;
                }

                Path mapPath = maps.get(fileIndex);
                String[] split = mapPath.toFile().getName().split("-");
                String stfId = split[0];

                String strategy = split[1];
                if(!strategy.equals("interacting")) {
                    continue;
                }

                // create directory
                stringJoiner.add("mkdir /tmp/reconstruction/");
                stringJoiner.add("mkdir /tmp/reconstruction/" + binIndex + "/");

                // create confold call
                stringJoiner.add("/home/bittrich/programs/confold_v1.0/confold.pl " +
                        "-rrtype ca " +
                        "-seq /home/bittrich/git/phd_sb_repo/data/start2fold/reconstruction/fasta/" + stfId + ".fasta " +
                        "-rr " + mapPath.toFile().getAbsolutePath() + " " +
                        "-ss /home/bittrich/git/phd_sb_repo/data/start2fold/reconstruction/sse/" + stfId + ".sse " +
                        "-o /tmp/reconstruction/" + binIndex + "/");

                String percentage = split[2];
                String samplingId = split[3].split("\\.")[0];
                // align to native structure
                // stfId, strategy, percentage, samplingId, confoldRunId
                stringJoiner.add("echo \"" + stfId + "," + strategy + "," + percentage + "," + samplingId + ",1\" >> confold-" + binIndex + "-interacting.out");
                stringJoiner.add("tmalign /tmp/reconstruction/" + binIndex + "/stage1/" + stfId + "_model1.pdb /home/bittrich/git/phd_sb_repo/data/start2fold/pdb/" + stfId + ".pdb | grep RMSD >> confold-" + binIndex + "-interacting.out");
                stringJoiner.add("echo \"" + stfId + "," + strategy + "," + percentage + "," + samplingId + ",2\" >> confold-" + binIndex + "-interacting.out");
                stringJoiner.add("tmalign /tmp/reconstruction/" + binIndex + "/stage1/" + stfId + "_model2.pdb /home/bittrich/git/phd_sb_repo/data/start2fold/pdb/" + stfId + ".pdb | grep RMSD >> confold-" + binIndex + "-interacting.out");
                stringJoiner.add("echo \"" + stfId + "," + strategy + "," + percentage + "," + samplingId + ",3\" >> confold-" + binIndex + "-interacting.out");
                stringJoiner.add("tmalign /tmp/reconstruction/" + binIndex + "/stage1/" + stfId + "_model3.pdb /home/bittrich/git/phd_sb_repo/data/start2fold/pdb/" + stfId + ".pdb | grep RMSD >> confold-" + binIndex + "-interacting.out");
                stringJoiner.add("echo \"" + stfId + "," + strategy + "," + percentage + "," + samplingId + ",4\" >> confold-" + binIndex + "-interacting.out");
                stringJoiner.add("tmalign /tmp/reconstruction/" + binIndex + "/stage1/" + stfId + "_model4.pdb /home/bittrich/git/phd_sb_repo/data/start2fold/pdb/" + stfId + ".pdb | grep RMSD >> confold-" + binIndex + "-interacting.out");
                stringJoiner.add("echo \"" + stfId + "," + strategy + "," + percentage + "," + samplingId + ",5\" >> confold-" + binIndex + "-interacting.out");
                stringJoiner.add("tmalign /tmp/reconstruction/" + binIndex + "/stage1/" + stfId + "_model5.pdb /home/bittrich/git/phd_sb_repo/data/start2fold/pdb/" + stfId + ".pdb | grep RMSD >> confold-" + binIndex + "-interacting.out");

                // cleanup
                stringJoiner.add("rm -r /tmp/reconstruction/" + binIndex + "/");
            }
            Start2FoldConstants.write(SCRIPT_DIRECTORY.resolve("confold-" + binIndex + "-interacting.sh"),
                    stringJoiner.toString());
        }

        for (int binIndex = 0; binIndex < bins; binIndex++) {
            StringJoiner stringJoiner = new StringJoiner(System.lineSeparator());
            for (int fileIndex = 0; fileIndex < maps.size(); fileIndex++) {
                if(fileIndex % bins != binIndex) {
                    continue;
                }

                Path mapPath = maps.get(fileIndex);
                String[] split = mapPath.toFile().getName().split("-");
                String stfId = split[0];

                String strategy = split[1];
                if(!strategy.equals("interacting2")) {
                    continue;
                }

                // create directory
                stringJoiner.add("mkdir /tmp/reconstruction/");
                stringJoiner.add("mkdir /tmp/reconstruction/" + binIndex + "/");

                // create confold call
                stringJoiner.add("/home/bittrich/programs/confold_v1.0/confold.pl " +
                        "-rrtype ca " +
                        "-seq /home/bittrich/git/phd_sb_repo/data/start2fold/reconstruction/fasta/" + stfId + ".fasta " +
                        "-rr " + mapPath.toFile().getAbsolutePath() + " " +
                        "-ss /home/bittrich/git/phd_sb_repo/data/start2fold/reconstruction/sse/" + stfId + ".sse " +
                        "-o /tmp/reconstruction/" + binIndex + "/");

                String percentage = split[2];
                String samplingId = split[3].split("\\.")[0];
                // align to native structure
                // stfId, strategy, percentage, samplingId, confoldRunId
                stringJoiner.add("echo \"" + stfId + "," + strategy + "," + percentage + "," + samplingId + ",1\" >> confold-" + binIndex + "-interacting2.out");
                stringJoiner.add("tmalign /tmp/reconstruction/" + binIndex + "/stage1/" + stfId + "_model1.pdb /home/bittrich/git/phd_sb_repo/data/start2fold/pdb/" + stfId + ".pdb | grep RMSD >> confold-" + binIndex + "-interacting2.out");
                stringJoiner.add("echo \"" + stfId + "," + strategy + "," + percentage + "," + samplingId + ",2\" >> confold-" + binIndex + "-interacting2.out");
                stringJoiner.add("tmalign /tmp/reconstruction/" + binIndex + "/stage1/" + stfId + "_model2.pdb /home/bittrich/git/phd_sb_repo/data/start2fold/pdb/" + stfId + ".pdb | grep RMSD >> confold-" + binIndex + "-interacting2.out");
                stringJoiner.add("echo \"" + stfId + "," + strategy + "," + percentage + "," + samplingId + ",3\" >> confold-" + binIndex + "-interacting2.out");
                stringJoiner.add("tmalign /tmp/reconstruction/" + binIndex + "/stage1/" + stfId + "_model3.pdb /home/bittrich/git/phd_sb_repo/data/start2fold/pdb/" + stfId + ".pdb | grep RMSD >> confold-" + binIndex + "-interacting2.out");
                stringJoiner.add("echo \"" + stfId + "," + strategy + "," + percentage + "," + samplingId + ",4\" >> confold-" + binIndex + "-interacting2.out");
                stringJoiner.add("tmalign /tmp/reconstruction/" + binIndex + "/stage1/" + stfId + "_model4.pdb /home/bittrich/git/phd_sb_repo/data/start2fold/pdb/" + stfId + ".pdb | grep RMSD >> confold-" + binIndex + "-interacting2.out");
                stringJoiner.add("echo \"" + stfId + "," + strategy + "," + percentage + "," + samplingId + ",5\" >> confold-" + binIndex + "-interacting2.out");
                stringJoiner.add("tmalign /tmp/reconstruction/" + binIndex + "/stage1/" + stfId + "_model5.pdb /home/bittrich/git/phd_sb_repo/data/start2fold/pdb/" + stfId + ".pdb | grep RMSD >> confold-" + binIndex + "-interacting2.out");

                // cleanup
                stringJoiner.add("rm -r /tmp/reconstruction/" + binIndex + "/");
            }
            Start2FoldConstants.write(SCRIPT_DIRECTORY.resolve("confold-" + binIndex + "-interacting2.sh"),
                    stringJoiner.toString());
        }
    }
}
