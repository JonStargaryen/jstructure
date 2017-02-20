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
        ObjectMapper mapper = new ObjectMapper();
        // ensure transient fields are not serialized
        mapper.configure(MapperFeature.PROPAGATE_TRANSIENT_MARKER, true);
        // TODO for development: prettify json
        mapper.enable(SerializationFeature.INDENT_OUTPUT);
        return mapper;
    }
}
