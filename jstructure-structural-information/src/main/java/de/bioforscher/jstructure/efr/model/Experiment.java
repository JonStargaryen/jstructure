package de.bioforscher.jstructure.efr.model;

import org.jsoup.nodes.Element;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class Experiment {
    private final int experimentId;
    private final Method method;
    private final ProtectionLevel protectionLevel;
    private final String sequence;
    private final List<Residue> residues;

    public Experiment(int experimentId,
                      Method method,
                      ProtectionLevel protectionLevel,
                      String sequence,
                      List<Residue> residues) {
        this.experimentId = experimentId;
        this.method = method;
        this.protectionLevel = protectionLevel;
        this.sequence = sequence;
        this.residues = residues;
    }

    public static Experiment parse(Element experiment) {
        return new Experiment(determineExperimentId(experiment),
                determineMethod(experiment),
                determineProtectionLevel(experiment),
                determineSequence(experiment),
                determineResidues(experiment));
    }

    private static int determineExperimentId(Element experiment) {
        return Integer.valueOf(experiment.attr("id"));
    }

    private static List<Residue> determineResidues(Element experiment) {
        return experiment.getElementsByTag("residue")
                .stream()
                .map(element -> new Residue(Integer.valueOf(element.attr("index")), element.attr("code")))
                .collect(Collectors.toList());
    }

    private static String determineSequence(Element experiment) {
        return experiment.getElementsByTag("sequence").first().text();
    }

    private static ProtectionLevel determineProtectionLevel(Element experiment) {
        String level = experiment.getElementsByTag("protection").first().attr("protection_level");
        return Stream.of(ProtectionLevel.values())
                .filter(protectionLevel -> protectionLevel.name().equalsIgnoreCase(level))
                .findFirst()
                .get();
    }

    private static Method determineMethod(Element experiment) {
        String type = experiment.getElementsByTag("method").first().attr("type");
        return type.equals("folding") ? Method.FOLDING : Method.STABILITY;
    }

    public int getExperimentId() {
        return experimentId;
    }

    public Method getMethod() {
        return method;
    }

    public ProtectionLevel getProtectionLevel() {
        return protectionLevel;
    }

    public String getSequence() {
        return sequence;
    }

    public List<Residue> getResidues() {
        return residues;
    }

    public static class Residue {
        private final int index;
        private final String code;

        public Residue(int index, String code) {
            this.index = index;
            this.code = code;
        }

        public int getIndex() {
            return index;
        }

        public String getCode() {
            return code;
        }
    }
}