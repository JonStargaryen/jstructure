package de.bioforscher.explorer;

import com.fasterxml.jackson.databind.MapperFeature;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.SerializationFeature;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;

@Configuration
public class JacksonConfiguration {
    @Bean
    public ObjectMapper objectMapper() {
        return new ObjectMapper()
        // ensure transient fields are not serialized
        .configure(MapperFeature.PROPAGATE_TRANSIENT_MARKER, true)
        // do not infer fields by getters
//        .configure(MapperFeature.AUTO_DETECT_GETTERS, false)
//        .configure(MapperFeature.AUTO_DETECT_IS_GETTERS, false)
        // TODO for development: prettify json
        .enable(SerializationFeature.INDENT_OUTPUT);
    }
}
