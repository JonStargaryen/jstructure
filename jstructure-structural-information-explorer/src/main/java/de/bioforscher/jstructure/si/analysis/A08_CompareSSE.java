package de.bioforscher.jstructure.si.analysis;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class A08_CompareSSE {
    public static void main(String[] args) throws IOException {
        Map<String, List<Double>> originalSR = parseSR("/Users/sebastian/1hrc-original.out");
        Map<String, List<Double>> originalRMSD = parseRMSD("/Users/sebastian/1hrc-original.out");
        Map<String, List<Double>> nosseSR = parseSR("/Users/sebastian/1hrc-no-sse-0.6.out");
        Map<String, List<Double>> nosseRMSD = parseRMSD("/Users/sebastian/1hrc-no-sse-0.6.out");

        String output = originalSR.keySet()
                .stream()
                .flatMap(key -> {
                    if (originalSR.get(key) == null || nosseSR.get(key) == null) {
                        return Stream.empty();
                    }

                    String[] split = key.replace("(", "")
                            .replace(")", "")
                            .split(", ");

//                    List<Double> sormsd = originalRMSD.get(key);
//                    List<Double> snrmsd = nosseRMSD.get(key);
//                    List<Double> sosr = originalSR.get(key);
//                    List<Double> snsr = nosseSR.get(key);
//
//                    Collections.sort(sormsd);
//                    Collections.sort(snrmsd);
//                    Collections.sort(sosr);
//                    Collections.sort(snsr);
//
//                    System.out.println(split[0] + " " + split[1] + " " + sormsd.size() + " " + snrmsd.size() + " " + sosr.size() + " " + snsr.size());
//
//                    return IntStream.range(0, Math.min(sormsd.size(), snrmsd.size()) - 1)
//                            .mapToObj(i -> split[0] + "," + split[1] + "," + sormsd.get(i) + "," + snrmsd.get(i) + "," + sosr.get(i) + "," + snsr.get(i));
                    return Stream.of(split[0] + "," +
                            split[1] + "," +
                            originalRMSD.get(key).stream().mapToDouble(d -> d).average().orElse(0) + "," +
                            nosseRMSD.get(key).stream().mapToDouble(d -> d).average().orElse(0) + "," +
                            originalSR.get(key).stream().mapToDouble(d -> d).average().orElse(0) + "," +
                            nosseSR.get(key).stream().mapToDouble(d -> d).average().orElse(0));
                })
                .collect(Collectors.joining("\n", "res1,res2,rmsdOriginal,rmsdNosse,srOriginal,srNosse\n", ""));

        System.out.println(output);

        Files.write(Paths.get("/Users/sebastian/sse-vs-nosse-0.6.csv"), output.getBytes());
    }

    private static Map<String, List<Double>> parseSR(String path) throws IOException {
        Map<String, List<Double>> map = new HashMap<>();
        Files.lines(Paths.get(path))
                .forEach(line -> {
                    String[] split = line.split("\t");
                    String key = split[0];
                    List<Double> value = map.computeIfAbsent(key, k -> new ArrayList<>());
                    value.add(Double.parseDouble(split[8]));
                });
        return map;
    }

    private static Map<String, List<Double>> parseRMSD(String path) throws IOException {
        Map<String, List<Double>> map = new HashMap<>();
        Files.lines(Paths.get(path))
                .forEach(line -> {
                    String[] split = line.split("\t");
                    String key = split[0];
                    List<Double> value = map.computeIfAbsent(key, k -> new ArrayList<>());
                    value.add(Double.parseDouble(split[5]));
                });
        return map;
    }
}
