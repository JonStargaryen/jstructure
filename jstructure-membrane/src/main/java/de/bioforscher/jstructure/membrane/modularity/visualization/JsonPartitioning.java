package de.bioforscher.jstructure.membrane.modularity.visualization;

import java.util.List;

public class JsonPartitioning {
    private final String identifier;
    private final List<JsonModule> modules;

    public JsonPartitioning(String identifier, List<JsonModule> modules) {
        this.identifier = identifier;
        this.modules = modules;
    }

    public String getIdentifier() {
        return identifier;
    }

    public List<JsonModule> getModules() {
        return modules;
    }
}
