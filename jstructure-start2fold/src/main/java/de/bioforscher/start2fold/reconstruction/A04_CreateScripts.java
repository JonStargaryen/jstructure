package de.bioforscher.start2fold.reconstruction;

import de.bioforscher.start2fold.Start2FoldConstants;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

public class A04_CreateScripts {
    private static final Path SCRIPT_DIRECTORY = Paths.get("/home/bittrich/git/phd_sb_repo/data/start2fold/reconstruction/");

    public static void main(String[] args) throws IOException {
        List<Path> maps = Files.list(Paths.get("/home/bittrich/git/phd_sb_repo/data/start2fold/reconstruction/maps/"))
                .limit(8)
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

                // create directory
                stringJoiner.add("mkdir /tmp/reconstruction/" + binIndex + "/");

                // create confold call
                stringJoiner.add("/home/bittrich/programs/confold_v1.0/confold.pl " +
                        "-rrtype ca " +
                        "-seq /home/bittrich/git/phd_sb_repo/data/start2fold/reconstruction/fasta/" + stfId + ".fasta " +
                        "-rr " + mapPath.toFile().getAbsolutePath() + " " +
                        "-ss /home/bittrich/git/phd_sb_repo/data/start2fold/reconstruction/sse/" + stfId + ".sse " +
                        "-o /tmp/reconstruction/" + binIndex + "/");

                // align to native structure

            }
            Start2FoldConstants.write(SCRIPT_DIRECTORY.resolve("confold-" + binIndex + ".sh"),
                    stringJoiner.toString());
        }
    }
}
