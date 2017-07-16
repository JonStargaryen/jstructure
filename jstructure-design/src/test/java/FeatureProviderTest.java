import de.bioforscher.jstructure.feature.energyprofile.EnergyProfile;
import de.bioforscher.jstructure.feature.sse.GenericSecondaryStructure;
import de.bioforscher.jstructure.model.feature.AbstractFeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureContainerEntry;
import de.bioforscher.jstructure.model.feature.FeatureProvider;
import de.bioforscher.jstructure.model.feature.FeatureProviderRegistry;
import de.bioforscher.jstructure.model.structure.Structure;
import org.junit.Test;

/**
 * Test for the {@link FeatureProviderRegistry} from a different module.
 * Created by bittrich on 7/5/17.
 */
public class FeatureProviderTest {
    @Test
    public void shouldResolveProvider() {
        FeatureProviderRegistry.getRegisteredFeatureProviders()
                .forEach(System.out::println);
        FeatureProviderRegistry.resolve(GenericSecondaryStructure.class);
        FeatureProviderRegistry.resolve(EnergyProfile.class);
    }

    @Test
    public void shouldRegisterAndResolveNewProvider() {
        //TODO other modules cannot register additional feature providers - using forceRegister as workaround
        FeatureProviderRegistry.register(AdditionalFeatureProvider.class);
        FeatureProviderRegistry.resolve(AdditionalFeatureEntry.class);
    }

    public static class AdditionalFeatureEntry extends FeatureContainerEntry {
        public AdditionalFeatureEntry(AbstractFeatureProvider featureProvider) {
            super(featureProvider);
        }
    }

    @FeatureProvider(provides = AdditionalFeatureEntry.class)
    public static class AdditionalFeatureProvider extends AbstractFeatureProvider {
        @Override
        protected void processInternally(Structure protein) {

        }
    }
}
