package de.bioforscher.jstructure.mmm.impl;

import de.bioforscher.jstructure.StandardFormat;
import de.bioforscher.jstructure.mmm.MacromolecularMinerBridge;
import de.bioforscher.jstructure.mmm.StructureConservationProfile;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.model.structure.StructureParser;
import de.bioforscher.testutil.TestUtils;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.concurrent.ExecutionException;
import java.util.stream.Stream;

/**
 * Tests for the macromolecular miner integration.
 * Created by bittrich on 7/13/17.
 */
public class MacromolecularMinerBridgeImplTest {
    private static final Logger logger = LoggerFactory.getLogger(MacromolecularMinerBridgeImplTest.class);
    private MacromolecularMinerBridge mmm;
    private Path in;
    private Path out;

    @Before
    public void setup() throws IOException {
        mmm = new MacromolecularMinerBridgeImpl();
        in = Files.createTempDirectory("mmm-in");
        out = Files.createTempDirectory("mmm-out");
        // fail-fast check for correct mmm integration
        mmm.getConservationProfileConfiguration(in, out);
        String ensemble = "1cyo_A, 1ehb_A, 1es1_A, 1f03_A, 1f04_A, 1hko_A, 1i5u_A, 1j0q_A, 1lqx_A, 1lr6_A, 1m20_A, 1m2i_A, 1m2m_A, 1m59_A, 1nx7_A, 1sh4_A, 1u9m_A, 1u9m_B, 1u9m_C, 1u9m_D, 1u9m_E, 1u9m_F, 1u9u_A, 3ozz_B, 4hin_A, 4hin_B, 4hin_C, 4hin_D, 3x32_A, 3x33_A, 3x34_A, 3x35_A, 1do9_A, 2m33_A, 2i96_A, 1aqa_A, 1aw3_A, 1axx_A, 1b5a_A, 1b5b_A, 1bfx_A, 1blv_A, 1i87_A, 1i8c_A, 1ib7_A, 1iet_A, 1ieu_A, 1jex_A, 1mny_A, 2axx_A, 1awp_A, 1awp_B, 1b5m_A, 1eue_A, 1eue_B, 1icc_A, 1icc_B, 1icc_C, 1icc_D, 1lj0_A, 1lj0_B, 1lj0_C, 1lj0_D, 2i89_A, 2i89_B, 2i89_C, 2i89_D, 3mus_A, 3mus_B, 4hil_A, 4hil_B, 3ner_A, 3ner_B, 2ibj_A, 1cne_A, 1cnf_A, 2cnd_A, 4zr0_A, 4zr0_B, 4zr1_A, 4zr1_B, 3lf5_A, 3lf5_B, 1cxy_A, 1fcb_A, 1fcb_B, 1kbi_A, 1kbi_B, 1kbj_A, 1kbj_B, 1lco_A, 1lco_B, 1ldc_A, 1ldc_B, 1ltd_A, 1ltd_B, 1qcw_A, 1qcw_B, 1sze_A, 1sze_B, 1szf_A, 1szf_B, 1szg_A, 1szg_B, 2oz0_A, 2oz0_B, 3ks0_A, 3ks0_B, 1mj4_A, 2bih_A, 2bii_A, 2bii_B, 1sox_A, 1sox_B, 2a99_A, 2a9a_A, 2a9a_B, 2a9b_A, 2a9c_A, 2a9c_B, 2a9d_A, 2a9d_B, 3hbg_A, 3hbp_A, 3hbq_A, 3r18_A, 3r19_A";
        // download files of ensemble and move to temp dir
        Stream.of(ensemble.split(", "))
                .map(string -> string.substring(0, 4))
                .distinct()
                .peek(pdbId -> logger.info("downloading {} of ensemble", pdbId))
                .forEach(pdbId -> TestUtils.downloadPdbId(in, pdbId));
    }

    @Test
    @Ignore
    public void shouldExecuteNormalMiningTask() throws ExecutionException, InterruptedException {
        mmm.submitStandardJob(in, out).get().getItemsetMiner();
    }

    @Test
    @Ignore
    public void shouldMineConservationProfile() throws ExecutionException, InterruptedException {
        Structure referenceProtein = StructureParser.source(TestUtils.getProteinInputStream(TestUtils.SupportedProtein.PDB_1CYO))
                .minimalParsing(true)
                .parse();
        mmm.getConservationProfile(in, referenceProtein).get();
        referenceProtein.aminoAcids()
                .map(aminoAcid -> aminoAcid.getIdentifier() + ": " + StandardFormat.format(aminoAcid.getFeatureContainer().getFeature(StructureConservationProfile.class).getValue()))
                .forEach(logger::info);
    }
}