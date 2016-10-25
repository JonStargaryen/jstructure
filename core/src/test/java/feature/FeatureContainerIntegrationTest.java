package feature;

import de.bioforscher.jstructure.feature.FeatureContainer;
import org.junit.Assert;
import org.junit.Test;

import java.util.HashMap;
import java.util.Map;

import static feature.FeatureContainerIntegrationTest.FeatureNames.*;

/**
 * Checks whether feature container capabilities are working.
 * Created by S on 25.10.2016.
 */
public class FeatureContainerIntegrationTest {
    private FeatureContainer featureContainer = new FeatureContainer() {
        private Map<Enum, Object> featureMap = new HashMap<>();

        @Override
        public Map<Enum, Object> getFeatureMap() {
            return featureMap;
        }
    };

    private static final int INT_VALUE_1 = 5;
    private static final int INT_VALUE_2 = 2;
    private static final double DOUBLE_VALUE = Math.PI;

    enum FeatureNames {
        PUT_INT,
        PUT_DOUBLE,
        PUT_OBJECT
    }

    @Test
    public void shouldPutElementIntoContainer() {
        featureContainer.setFeature(PUT_INT, INT_VALUE_1);
        featureContainer.setFeature(PUT_DOUBLE, DOUBLE_VALUE);
        featureContainer.setFeature(PUT_OBJECT, this);
        Assert.assertEquals(INT_VALUE_1, featureContainer.getIntFeature(PUT_INT));
        Assert.assertEquals(DOUBLE_VALUE, featureContainer.getDoubleFeature(PUT_DOUBLE), 0.0);
        Assert.assertEquals(this, featureContainer.getFeature(getClass(), PUT_OBJECT));
    }

    @Test
    public void shouldReplaceElementInContainer() {
        featureContainer.setFeature(PUT_INT, INT_VALUE_1);
        Assert.assertEquals(INT_VALUE_1, featureContainer.getIntFeature(PUT_INT));
        featureContainer.setFeature(PUT_INT, INT_VALUE_2);
        Assert.assertEquals(INT_VALUE_2, featureContainer.getIntFeature(PUT_INT));
    }
}
