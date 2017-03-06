package de.bioforscher.explorer;

/**
 * Allows to 'mock' the persistence layer when none is available.
 * TODO find better pattern to switch
 * Created by S on 04.03.2017.
 */
//@RestController
//@RequestMapping(value = "/api/proteins/", method = RequestMethod.GET)
public class MockController {
//    @RequestMapping(value = "/all", method = RequestMethod.GET)
//    public List<String> getAllProteins() throws IOException {
//        return Stream.of("5A63").collect(Collectors.toList());
//    }
//
//    @RequestMapping(value = "/alpha_nr", method = RequestMethod.GET)
//    public List<String> getNonRedundantAlphaHelicalProteins() throws IOException {
//        return Stream.of("5A63").collect(Collectors.toList());
//    }
//
//    @RequestMapping(value = "/pdbid/{pdbid}", method = RequestMethod.GET)
//    public String getProteinByID(@PathVariable("pdbid") String pdbid) throws IOException {
//        return Files.lines(getResource("mock.json")).collect(Collectors.joining(System.lineSeparator()));
//    }
//
//    private Path getResource(String filename) {
//        // some a bit hacky way to ensure correct paths on windows (as some / will be added as prefix)
//        return Paths.get(Objects.requireNonNull(Thread.currentThread().getContextClassLoader().getResource(filename)).getPath().replaceFirst("^/(.:/)", "$1"));
//    }
}
