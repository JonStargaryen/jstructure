package de.bioforscher.jstructure.mutation.impl;

import de.bioforscher.jstructure.model.feature.*;
import de.bioforscher.jstructure.model.identifier.ChainIdentifier;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.aminoacid.AminoAcid;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Query the ConSurf for scores.
 * Created by bittrich on 7/19/17.
 */
@FeatureProvider(provides = ConSurfDBQuery.ConSurfScore.class)
public class ConSurfDBQuery extends AbstractFeatureProvider {
    private static final Logger logger = LoggerFactory.getLogger(ConSurfDBQuery.class);
    private static final String CONSURF_QUERY_URL = "http://bental.tau.ac.il/new_ConSurfDB/main_output.php?pdb_ID=%s&view_chain=files_%s";

    protected void processInternally(Structure protein) {
        protein.chainsWithAminoAcids().forEach(this::processInternally);
    }

    private void processInternally(Chain chain) {
        ChainIdentifier chainIdentifier = chain.getChainIdentifier();
        try {
            // use chain-specific entry to navigate to representative
            Document document = Jsoup.connect(String.format(CONSURF_QUERY_URL, chainIdentifier.getProteinIdentifier().getPdbId(), chainIdentifier.getChainId())).get();

            String fileLocation = document.getElementsByTag("form").stream()
                    .map(element -> element.attr("action"))
                    .filter(action -> action.contains("consurf.grades"))
                    .findFirst()
                    .get();

            // parse normalized scores
            List<Double> scores = parseGradesFiles(new URL(fileLocation).openStream());
            List<AminoAcid> aminoAcids = chain.aminoAcids().collect(Collectors.toList());

            for(int i = 0; i < aminoAcids.size(); i++) {
                aminoAcids.get(i).getFeatureContainer().addFeature(new ConSurfScore(this, scores.get(i)));
            }
        } catch (Exception e) {
            logger.warn("failed to fetch ConSurf data for chain {}",
                    chainIdentifier.getFullName(),
                    e);
            throw new ComputationException(e);
        }
    }

    public List<Double> parseGradesFiles(InputStream inputStream) {
        return new BufferedReader(new InputStreamReader(inputStream))
                .lines()
                .map(String::trim)
                // this criteria breaks easily
                .filter(line -> line.contains(",") && line.contains("/"))
                .map(line -> line.split("\\s+"))
                .map(split -> split[3])
                .map(Double::valueOf)
                .collect(Collectors.toList());
    }

    public static class ConSurfScore extends FeatureContainerEntry implements SingleValueFeatureContainerEntry<Double> {
        private final double score;

        ConSurfScore(AbstractFeatureProvider featureProvider, double score) {
            super(featureProvider);
            this.score = score;
        }

        public double getScore() {
            return score;
        }

        @Override
        public Double getValue() {
            return score;
        }
    }
}
