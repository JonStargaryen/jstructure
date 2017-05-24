package de.bioforscher.jstructure.model.structure.prototype;

import java.io.*;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Map;
import java.util.Objects;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Parses CIF files locally stored for standard-cases or fetches them from the PDB.
 * Created by bittrich on 5/24/17.
 */
class PrototypeParser {
    /**
     * The URL which can be used to fetch ligand information.
     */
    private static final String CIF_FETCH_URL = "https://files.rcsb.org/ligands/view/%s.cif";
    private static final PrototypeParser INSTANCE = new PrototypeParser();
    /**
     * Collection of already parsed group information.
     */
    private final Map<String, GroupPrototype> prototypes;

    private PrototypeParser() {
        try {
            this.prototypes = Files.list(Paths.get(getResourceAsFilepath("prototype/")))
                    .map(this::createPrototype)
                    .collect(Collectors.toConcurrentMap(GroupPrototype::getId,
                            Function.identity()));
        } catch (IOException e) {
            throw new UncheckedIOException("could not find local collection of CIF files, failed to initialize parsers", e);
        }
    }

    private GroupPrototype createPrototype(Path path) {
        try {
            return createPrototype(Files.newInputStream(path));
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    private GroupPrototype createPrototype(InputStream inputStream) throws IOException {
        try(InputStreamReader inputStreamReader = new InputStreamReader(inputStream)) {
            try (BufferedReader bufferedReader = new BufferedReader(inputStreamReader)) {

            }
        }
    }

    private GroupPrototype createPrototype(String id) {
        //TODO continue here
    }

    private String getResourceAsFilepath(String filename) {
        ClassLoader ccl = Thread.currentThread().getContextClassLoader();
        URL resource = ccl.getResource(filename);
        Objects.requireNonNull(resource);
        // some a bit hacky way to ensure correct paths on windows (as some / will be added as prefix)
        return resource.getPath().replaceFirst("^/(.:/)", "$1");
    }

    public static PrototypeParser getInstance() {
        return INSTANCE;
    }

    public GroupPrototype getPrototype(String id) {
        return prototypes.computeIfAbsent(id, this::createPrototype);
    }
}
