package de.bioforscher.jstructure.mmm.impl;

import de.bioforscher.jstructure.mmm.MacromolecularMinerBridge;
import de.bioforscher.jstructure.mmm.StructureConservationProfile;
import de.bioforscher.testutil.TestUtils;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

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
        String ensemble = "1aey_A, 1aj3_A, 1bk2_A, 1cun_A, 1cun_B, 1cun_C, 1e6g_A, 1e6h_A, 1e7o_A, 1g2b_A, 1h8k_A, 1hd3_A, 1m8m_A, 1neg_A, 1pwt_A, 1qkw_A, 1qkx_A, 1shg_A, 1tuc_A, 1tud_A, 1u06_A, 1u4q_A, 1u4q_B, 1u5p_A, 1uue_A, 2cdt_A, 2f2v_A, 2f2w_A, 2f2x_A, 2jm8_A, 2jm9_A, 2jma_A, 2jmc_A, 2kr3_A, 2kxd_A, 2lj3_A, 2nuz_A, 2oaw_A, 2oaw_B, 2oaw_C, 2oaw_D, 2rmo_A, 2rot_A, 3i9q_A, 3m0p_A, 3m0q_A, 3m0r_A, 3m0s_A, 3m0t_A, 3m0u_A, 3ngp_A, 4f16_A, 4f17_A, 5ihi_A, 5ihk_A, 5ihn_A, 5m6s_A, 3thk_A, 3thk_B, 2fot_C, 3f31_A, 3f31_B, 3fb2_A, 3fb2_B, 5fw9_A, 5fwb_A, 5fwc_A, 1owa_A, 3lbx_A, 5j4o_A, 2spc_A, 2spc_B, 1dro_A, 1btn_A, 1mph_A, 1aa2_A, 1bkr_A, 3edv_A, 3edv_B, 1s35_A, 3edu_A, 3f57_A, 3f57_B, 3kbt_A, 3kbt_B, 3kbu_A, 3kbu_B, 3lbx_B, 1wjm_A, 1wyq_A, 1sjj_A, 1sjj_B, 1h8b_A, 1hci_A, 1hci_B, 1quu_A, 4d1e_A, 5a36_A, 5a36_B, 5a37_A, 5a37_B, 5a38_A, 5a38_B, 5a4b_A, 5a4b_B, 2eyi_A, 2eyn_A, 2n8y_A, 2n8z_A";
        // download files of ensemble and move to temp dir
        Stream.of(ensemble.split(", "))
                .map(string -> string.substring(0, 4))
                .limit(10)
                .distinct()
                .forEach(pdbId -> TestUtils.downloadPdbId(in, pdbId));
    }

    @Test
    @Ignore
    public void shouldExecuteNormalMiningTask() throws IOException, ExecutionException, InterruptedException {
        mmm.submitStandardJob(in, out).get().getItemsetMiner();
    }

    @Test
    @Ignore
    public void shouldMineConservationProfile() throws IOException, ExecutionException, InterruptedException {
        StructureConservationProfile structureConservationProfile = mmm.getConservationProfile(in).get();
    }
}