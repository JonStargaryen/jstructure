package de.bioforscher.jstructure.contacts.collect.reconstruction.sathyapriya;

import de.bioforscher.jstructure.contacts.ContactsConstants;
import de.bioforscher.jstructure.contacts.collect.reconstruction.ContactMapCreator;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.StringJoiner;
import java.util.stream.Collectors;

public class A04_CreateConfoldScript {
    private static final Path DIRECTORY = ContactsConstants.RECONSTRUCTION_DIRECTORY;

    public static void main(String[] args) {
        String fullJobOutput = ContactsConstants.lines(ContactsConstants.RECONSTRUCTION_DIRECTORY.resolve("ids.list"))
                .map(A04_CreateConfoldScript::handlePdbIdForFullJob)
                .collect(Collectors.joining(System.lineSeparator()));

        ContactsConstants.write(DIRECTORY.resolve("full-reconstruction.sh"), fullJobOutput);

        for(ContactMapCreator.Sampling sampling : ContactMapCreator.Sampling.values()) {
            String sampledJobOutput = ContactsConstants.lines(ContactsConstants.RECONSTRUCTION_DIRECTORY.resolve("ids.list"))
                    .map(pdbId -> handlePdbIdForSampledJob(pdbId, sampling.getFraction()))
                    .collect(Collectors.joining(System.lineSeparator()));

            ContactsConstants.write(DIRECTORY.resolve("sampled-" + sampling.getFraction() + "-reconstruction.sh"), sampledJobOutput);
        }
    }

    private static String handlePdbIdForSampledJob(String pdbId, int sampling) {
        StringJoiner stringJoiner = new StringJoiner(System.lineSeparator());
        // select 10 maps of each run
        for (int i = 1; i < 11; i++) {
            try {
                Files.createDirectories(Paths.get("/home/bittrich/git/phd_sb_repo/data/reconstruction/reconstructions/p" + sampling + "/" + pdbId + "-naive-" + i + "/"));
                Files.createDirectories(Paths.get("/home/bittrich/git/phd_sb_repo/data/reconstruction/reconstructions/p" + sampling + "/" + pdbId + "-plip-" + i + "/"));
                stringJoiner.add("/home/bittrich/programs/confold_v1.0/confold.pl -seq /home/bittrich/git/phd_sb_repo/data/reconstruction/maps/" + pdbId + "_A.fasta -rr /home/bittrich/git/phd_sb_repo/data/reconstruction/maps/p" + sampling + "/" + pdbId + "_A-naive-" + i + ".rr -o /home/bittrich/git/phd_sb_repo/data/reconstruction/reconstructions/p" + sampling + "/" + pdbId + "-naive-" + i + "/");
                stringJoiner.add("/home/bittrich/programs/confold_v1.0/confold.pl -seq /home/bittrich/git/phd_sb_repo/data/reconstruction/maps/" + pdbId + "_A.fasta -rr /home/bittrich/git/phd_sb_repo/data/reconstruction/maps/p" + sampling + "/" + pdbId + "_A-plip-" + i + ".rr -o /home/bittrich/git/phd_sb_repo/data/reconstruction/reconstructions/p" + sampling + "/" + pdbId + "-plip-" + i + "/");
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
        }
        return stringJoiner.toString();
    }

    private static String handlePdbIdForFullJob(String pdbId) {
        try {
            Files.createDirectories(Paths.get("/home/bittrich/git/phd_sb_repo/data/reconstruction/reconstructions/p100/" + pdbId + "-naive/"));
            Files.createDirectories(Paths.get("/home/bittrich/git/phd_sb_repo/data/reconstruction/reconstructions/p100/" + pdbId + "-plip/"));
            return "/home/bittrich/programs/confold_v1.0/confold.pl -seq /home/bittrich/git/phd_sb_repo/data/reconstruction/maps/" + pdbId + "_A.fasta -rr /home/bittrich/git/phd_sb_repo/data/reconstruction/maps/p100/" + pdbId + "_A-naive.rr -o /home/bittrich/git/phd_sb_repo/data/reconstruction/reconstructions/p100/" + pdbId + "-naive/" + System.lineSeparator() +
                    "/home/bittrich/programs/confold_v1.0/confold.pl -seq /home/bittrich/git/phd_sb_repo/data/reconstruction/maps/" + pdbId + "_A.fasta -rr /home/bittrich/git/phd_sb_repo/data/reconstruction/maps/p100/" + pdbId + "_A-plip.rr -o /home/bittrich/git/phd_sb_repo/data/reconstruction/reconstructions/p100/" + pdbId + "-plip/";
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }
}
