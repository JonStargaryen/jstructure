package aars;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

/**
 * Compose weblogo-able groups of strings for the backbone-brackets as well as arginine tweezers.
 *
 * http://weblogo.berkeley.edu/logo.cgi
 *
 * settings:
 * image format: PDF
 * position numbers: provide them
 * show sequence ends: false
 * show fine print: false
 * color scheme: custom
 * c, S, T - coil - black
 * I, G, H - helix - blue
 * B, E - strand - red
 *
 * Created by bittrich on 4/19/17.
 */
public class CreateSecondaryStructureWebLogos {
    public static void main(String[] args) throws IOException {
        handleFile("/home/bittrich/git/aars_analysis/data/geometry/argtweezer_geometry.tsv");
        handleFile("/home/bittrich/git/aars_analysis/data/geometry/bbbrackets_geometry.tsv");
    }

    private static void handleFile(String filepath) throws IOException {
        System.out.println(filepath);

        System.out.println("position 1:");
        Files.lines(Paths.get(filepath))
                .filter(line -> !line.startsWith("id"))
                .map(line -> line.split("\t"))
                .map(split -> split[5])
                .forEach(System.out::println);
        System.out.println();

        System.out.println("position 2:");
        Files.lines(Paths.get(filepath))
                .filter(line -> !line.startsWith("id"))
                .map(line -> line.split("\t"))
                .map(split -> split[6])
                .forEach(System.out::println);
        System.out.println();
    }
}
