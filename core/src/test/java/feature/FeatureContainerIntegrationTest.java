package feature;

import de.bioforscher.jstructure.model.feature.AbstractFeatureContainer;
import de.bioforscher.jstructure.model.feature.FeatureContainer;
import org.junit.Assert;
import org.junit.Test;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static feature.FeatureContainerIntegrationTest.FeatureNames.*;

/**
 * Checks whether feature container capabilities are working.
 * Created by S on 25.10.2016.
 */
public class FeatureContainerIntegrationTest {
    private FeatureContainer featureContainer = new AbstractFeatureContainer() {
        @Override
        public String getIdentifier() {
            return null;
        }

        @Override
        public void setIdentifier(String identifier) {

        }

        @Override
        public String composePDBRecord() {
            return null;
        }
    };

    private static final int INT_VALUE_1 = 5;
    private static final int INT_VALUE_2 = 2;
    private static final double DOUBLE_VALUE = Math.PI;

    enum FeatureNames {
        PUT_INT,
        PUT_DOUBLE,
        PUT_OBJECT,
        PUT_LIST
    }

    @Test
    public void shouldPutElementIntoContainer() {
        featureContainer.setFeature(PUT_INT, INT_VALUE_1);
        featureContainer.setFeature(PUT_DOUBLE, DOUBLE_VALUE);
        featureContainer.setFeature(PUT_OBJECT, this);
        Assert.assertEquals(INT_VALUE_1, featureContainer.getFeatureAsInt(PUT_INT));
        Assert.assertEquals(DOUBLE_VALUE, featureContainer.getFeatureAsDouble(PUT_DOUBLE), 0.0);
        Assert.assertEquals(this, featureContainer.getFeature(getClass(), PUT_OBJECT));
    }

    @Test
    public void shouldReplaceElementInContainer() {
        featureContainer.setFeature(PUT_INT, INT_VALUE_1);
        Assert.assertEquals(INT_VALUE_1, featureContainer.getFeatureAsInt(PUT_INT));
        featureContainer.setFeature(PUT_INT, INT_VALUE_2);
        Assert.assertEquals(INT_VALUE_2, featureContainer.getFeatureAsInt(PUT_INT));
    }

    @Test
    public void testBehaviourWithLists() {
        List<Double> values = Stream.of(1.0, 2.0, 3.0, 4.0, 5.0).collect(Collectors.toList());
        featureContainer.setFeature(PUT_LIST, values);

        List<Double> content = featureContainer.getFeatureAsList(Double.class, PUT_LIST);
        content.stream()
                .peek(System.out::println)
                .map(Double.class::isInstance)
                .forEach(Assert::assertTrue);
    }

    @Test(expected = ClassCastException.class)
    public void testBehaviourWithListsForErroneousCast_expectClassCastException() {
        List<Double> values = Stream.of(1.0, 2.0, 3.0, 4.0, 5.0).collect(Collectors.toList());
        featureContainer.setFeature(PUT_LIST, values);

        // wrong content type provided, will throw ClassCastException
        featureContainer.getFeatureAsList(Boolean.class, PUT_LIST);
    }
}
