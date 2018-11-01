package de.bioforscher.jstructure.feature.plip;

import de.bioforscher.jstructure.feature.plip.model.*;
import de.bioforscher.jstructure.model.feature.ComputationException;
import de.bioforscher.jstructure.model.structure.Atom;
import de.bioforscher.jstructure.model.structure.Chain;
import de.bioforscher.jstructure.model.structure.Group;
import de.bioforscher.jstructure.model.structure.Structure;
import de.bioforscher.jstructure.service.ExternalLocalService;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.nodes.Element;
import org.jsoup.select.Elements;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.net.URLConnection;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Annotates interactions between macromolecules and ligands as well as polymer interaction within macromolecules.
 *
 * Salentin,S. et al. PLIP: fully automated protein-ligand interaction profiler. "Nucl. Acids Res. (1 July 2015) 43
 * (W1): W443-W447. doi: 10.1093/nar/gkv315
 */
public class ProteinLigandInteractionProfiler extends ExternalLocalService {
    private static final Logger logger = LoggerFactory.getLogger(ProteinLigandInteractionProfiler.class);
    private static final ProteinLigandInteractionProfiler INSTANCE = new ProteinLigandInteractionProfiler();
    private static final String DEFAULT_SERVICE_LOCATION = "/home/bittrich/programs/pliptool/plip/plipcmd.py";
    private static final String BASE_URL = "https://biosciences.hs-mittweida.de/plip/api/";
    private String secret;

    private ProteinLigandInteractionProfiler() {
        super(DEFAULT_SERVICE_LOCATION);
        try {
            String line = new BufferedReader(new InputStreamReader(Thread.currentThread()
                    .getContextClassLoader()
                    .getResourceAsStream("plip_credentials.txt"))).readLine();
            secret = new String(Base64.getMimeEncoder().encode(line.getBytes()));
            logger.debug("PLIP accessed via {} - authentication provided", BASE_URL);
        } catch (IOException | NullPointerException e) {
            throw new IllegalStateException("no credentials provided to access 'biosciences.hs-mittweida.de/plip/'");
        }
    }

    public static ProteinLigandInteractionProfiler getInstance() {
        return INSTANCE;
    }

    public ProteinLigandInteractionProfiler(String serviceLocation) {
        super(serviceLocation);
    }

    public InteractionContainer annotate(Structure structure) {
        List<Interaction> mergedInteractions = annotateLigandInteractions(structure).getInteractions();
        mergedInteractions.addAll(annotatePolymerInteractions(structure).getInteractions());
        return new InteractionContainer(mergedInteractions);
    }

    public InteractionContainer fetchLigandInteractions(Structure structure) {
        try {
            String pdbId = structure.getProteinIdentifier().getPdbId();
            String boundary = Long.toHexString(System.currentTimeMillis());
            URLConnection connection = (new URL("https://biosciences.hs-mittweida.de/plip/api/provided/ligand/" + pdbId)).openConnection();
            connection.setDoOutput(true);
            connection.setRequestProperty("Content-Type", "multipart/form-data; boundary=" + boundary);
            connection.setRequestProperty("Authorization", "Basic " + secret);
            BufferedReader in = new BufferedReader(new InputStreamReader(connection.getInputStream()));
            StringBuilder response = new StringBuilder();

            String outputContent;
            while((outputContent = in.readLine()) != null) {
                response.append(outputContent);
            }

            in.close();

            Document document = Jsoup.parse(response.toString());
            return parseDocument(structure, document);
        } catch (Exception e) {
            //TODO this happens as no ligand interactions are precomputed at biosciences
            throw new ComputationException(e);
        }
    }

    public InteractionContainer annotateLigandInteractions(Structure structure) {
        try {
            Path outputDirectoryPath = createTemporaryOutputDirectory();
            String[] commandLineCall = composeLigandCommandLineCall(structure,
                    outputDirectoryPath);

            // execute command
            executeCommandLineCall(commandLineCall);

            Path outputPath = outputDirectoryPath.resolve("report.xml");
            Document document = Jsoup.parse(outputPath.toFile(), "UTF-8");
            return parseDocument(structure, document);
        } catch (Exception e) {
            throw new ComputationException(e);
        }
    }

    public InteractionContainer annotatePolymerInteractions(Structure structure) {
        //TODO slow
        //TODO ensure uniqueness of interactions
        return new InteractionContainer(structure.chainsWithAminoAcids()
                .map(this::annotatePolymerInteractions)
                .map(InteractionContainer::getInteractions)
                .flatMap(Collection::stream)
                .collect(Collectors.toList()));
    }

    private InteractionContainer annotatePolymerInteractions(Chain chain) {
        try {
            Path outputDirectoryPath = createTemporaryOutputDirectory();
            Structure structure = chain.getParentStructure();
            String[] commandLineCall = composePolymerCommandLineCall(structure,
                    outputDirectoryPath,
                    chain.getChainIdentifier().getChainId());

            // execute command
            executeCommandLineCall(commandLineCall);

            Path outputPath = outputDirectoryPath.resolve("report.xml");
            Document document = Jsoup.parse(outputPath.toFile(), "UTF-8");
            return parseDocument(structure, document);
        } catch (Exception e) {
            throw new ComputationException(e);
        }
    }

    private InteractionContainer parseDocument(Structure structure, Document document) {
        List<Atom> filteredAtoms = structure.atoms()
                .filter(atom -> atom.getElement().isHeavyAtom())
                .filter(atom -> !atom.hasAlternativeLocations())
                .collect(Collectors.toList());
        List<Interaction> interactions = new ArrayList<>();
        for(Element elementTag : document.getElementsByTag("interactions")) {
            parseHalogenBonds(structure, elementTag, interactions);
            parseHydrogenBonds(structure, elementTag, interactions);
            parseHydrophobicInteractions(structure, elementTag, interactions);
            parseMetalComplexes(structure, elementTag, interactions);
            parsePiCationInteractions(structure, elementTag, interactions);
            parsePiStackingInteractions(structure, elementTag, interactions);
            parseSaltBridges(structure, elementTag, interactions);
            parseWaterBridges(structure, elementTag, interactions);
//            parseHalogenBonds(filteredAtoms, elementTag, interactions);
//            parseHydrogenBonds(filteredAtoms, elementTag, interactions);
//            parseHydrophobicInteractions(filteredAtoms, elementTag, interactions);
//            parseMetalComplexes(filteredAtoms, elementTag, interactions);
//            parsePiCationInteractions(filteredAtoms, elementTag, interactions);
//            parsePiStackingInteractions(filteredAtoms, elementTag, interactions);
//            parseSaltBridges(filteredAtoms, elementTag, interactions);
//            parseWaterBridges(filteredAtoms, elementTag, interactions);
        }
        //TODO: check redundancy
        return new InteractionContainer(interactions);
    }

    private void parseHalogenBonds(Structure structure, Element interactionElement, List<Interaction> interactions) {
//        private void parseHalogenBonds(List<Atom> filteredAtoms, Element interactionElement, List<Interaction> interactions) {
        Elements halogenBonds = interactionElement.getElementsByTag("halogen_bond");
        halogenBonds.stream()
                .map(element -> parseHalogenBond(structure, element))
//                .map(element -> parseHalogenBond(filteredAtoms, element))
                .filter(Optional::isPresent)
                .map(Optional::get)
                .forEach(interactions::add);
    }

    private Optional<HalogenBond> parseHalogenBond(Structure structure, Element element) {
        Optional<Atom> optionalDonor = selectAtom(structure, element.getElementsByTag("don_idx").first().text());
        Optional<Atom> optionalAcceptor = selectAtom(structure, element.getElementsByTag("acc_idx").first().text());
//    private Optional<HalogenBond> parseHalogenBond(List<Atom> filteredAtoms, Element element) {
//        Optional<Atom> optionalDonor = selectAtom(filteredAtoms, element.getElementsByTag("don_idx").first().text());
//        Optional<Atom> optionalAcceptor = selectAtom(filteredAtoms, element.getElementsByTag("acc_idx").first().text());
        if(!(optionalDonor.isPresent() && optionalAcceptor.isPresent())) {
            return Optional.empty();
        }
        Atom donor = optionalDonor.get();
        Atom acceptor = optionalAcceptor.get();
        Group donorGroup = donor.getParentGroup();
        Group acceptorGroup = acceptor.getParentGroup();
        return Optional.of(new HalogenBond(donor,
                acceptor,
                donorGroup,
                acceptorGroup));
    }

    private void parseHydrogenBonds(Structure structure, Element interactionElement, List<Interaction> interactions) {
//    private void parseHydrogenBonds(List<Atom> filteredAtoms, Element interactionElement, List<Interaction> interactions) {
        Elements hydrogenBonds = interactionElement.getElementsByTag("hydrogen_bond");
        hydrogenBonds.stream()
                .map(element -> parseHydrogenBond(structure, element))
//                .map(element -> parseHydrogenBond(filteredAtoms, element))
                .filter(Optional::isPresent)
                .map(Optional::get)
                .forEach(interactions::add);
    }

    private Optional<HydrogenBond> parseHydrogenBond(Structure structure, Element element) {
//    private Optional<HydrogenBond> parseHydrogenBond(List<Atom> filteredAtoms, Element element) {
        Optional<Atom> optionalDonor = selectAtom(structure, element.getElementsByTag("donoridx").first().text());
        Optional<Atom> optionalAcceptor = selectAtom(structure, element.getElementsByTag("acceptoridx").first().text());
//        Optional<Atom> optionalDonor = selectAtom(filteredAtoms, element.getElementsByTag("donoridx").first().text());
//        Optional<Atom> optionalAcceptor = selectAtom(filteredAtoms, element.getElementsByTag("acceptoridx").first().text());
        if(!(optionalDonor.isPresent() && optionalAcceptor.isPresent())) {
            return Optional.empty();
        }
        Atom donor = optionalDonor.get();
        Atom acceptor = optionalAcceptor.get();
        Group donorGroup = donor.getParentGroup();
        Group acceptorGroup = acceptor.getParentGroup();
        return Optional.of(new HydrogenBond(donor,
                acceptor,
                donorGroup,
                acceptorGroup));
    }

    private void parseHydrophobicInteractions(Structure structure, Element interactionElement, List<Interaction> interactions) {
//    private void parseHydrophobicInteractions(List<Atom> filteredAtoms, Element interactionElement, List<Interaction> interactions) {
        Elements hydrophobicInteractions = interactionElement.getElementsByTag("hydrophobic_interaction");
        hydrophobicInteractions.stream()
                .map(element -> parseHydrophobicInteraction(structure, element))
//                .map(element -> parseHydrophobicInteraction(filteredAtoms, element))
                .filter(Optional::isPresent)
                .map(Optional::get)
                .forEach(interactions::add);
    }

    private Optional<HydrophobicInteraction> parseHydrophobicInteraction(Structure structure, Element element) {
        Optional<Atom> optionalAtom1 = selectAtom(structure, element.getElementsByTag("ligcarbonidx").first().text());
        Optional<Atom> optionalAtom2 = selectAtom(structure, element.getElementsByTag("protcarbonidx").first().text());
//    private Optional<HydrophobicInteraction> parseHydrophobicInteraction(List<Atom> filteredAtoms, Element element) {
//        Optional<Atom> optionalAtom1 = selectAtom(filteredAtoms, element.getElementsByTag("ligcarbonidx").first().text());
//        Optional<Atom> optionalAtom2 = selectAtom(filteredAtoms, element.getElementsByTag("protcarbonidx").first().text());
        if(!(optionalAtom1.isPresent() && optionalAtom2.isPresent())) {
            return Optional.empty();
        }
        Atom atom1 = optionalAtom1.get();
        Atom atom2 = optionalAtom2.get();
        Group group1 = atom1.getParentGroup();
        Group group2 = atom2.getParentGroup();
        return Optional.of(new HydrophobicInteraction(atom1,
                atom2,
                group1,
                group2));
    }

    private void parseMetalComplexes(Structure structure, Element interactionElement, List<Interaction> interactions) {
//    private void parseMetalComplexes(List<Atom> filteredAtoms, Element interactionElement, List<Interaction> interactions) {
        Elements metalComplexes = interactionElement.getElementsByTag("metal_complex");
        metalComplexes.stream()
                .map(element -> parseMetalComplex(structure, element))
//                .map(element -> parseMetalComplex(filteredAtoms, element))
                .filter(Optional::isPresent)
                .map(Optional::get)
                .forEach(interactions::add);
    }

    private Optional<MetalComplex> parseMetalComplex(Structure structure, Element element) {
        //TODO: potentially half interaction
        Optional<Atom> optionalAtom1 = selectAtom(structure, element.getElementsByTag("metal_idx").first().text());
        Optional<Atom> optionalAtom2 = selectAtom(structure, element.getElementsByTag("target_idx").first().text());
//    private Optional<MetalComplex> parseMetalComplex(List<Atom> filteredAtoms, Element element) {
//        //TODO: potentially half interaction
//        Optional<Atom> optionalAtom1 = selectAtom(filteredAtoms, element.getElementsByTag("metal_idx").first().text());
//        Optional<Atom> optionalAtom2 = selectAtom(filteredAtoms, element.getElementsByTag("target_idx").first().text());
        if(!(optionalAtom1.isPresent() && optionalAtom2.isPresent())) {
            return Optional.empty();
        }
        Atom atom1 = optionalAtom1.get();
        Atom atom2 = optionalAtom2.get();
        Group group1 = atom1.getParentGroup();
        Group group2 = atom2.getParentGroup();
        return Optional.of(new MetalComplex(atom1,
                atom2,
                group1,
                group2));
    }

    private void parsePiCationInteractions(Structure structure, Element interactionElement, List<Interaction> interactions) {
//    private void parsePiCationInteractions(List<Atom> filteredAtoms, Element interactionElement, List<Interaction> interactions) {
        Elements piCationInteractions = interactionElement.getElementsByTag("pi_cation_interaction");
        piCationInteractions.stream()
                .map(element -> parsePiCationInteraction(structure, element))
//                .map(element -> parsePiCationInteraction(filteredAtoms, element))
                .filter(Optional::isPresent)
                .map(Optional::get)
                .forEach(interactions::add);
    }

    private Optional<PiCationInteraction> parsePiCationInteraction(Structure structure, Element element) {
//    private Optional<PiCationInteraction> parsePiCationInteraction(List<Atom> filteredAtoms, Element element) {
//        Structure structure = filteredAtoms.get(0).getParentStructure();

        //TODO check if direction can be reversed: protein is pi system, ligand is cation
        //TODO: potentially half interaction
        Optional<Element> optionalCationCoordinateElement = element.getElementsByTag("protcoo").stream().findFirst();
        if(!optionalCationCoordinateElement.isPresent()) {
            return Optional.empty();
        }

        Element cationCoordinateElement = optionalCationCoordinateElement.get();
        double[] cationCoordinates = new double[] {
                Double.valueOf(cationCoordinateElement.getElementsByTag("x").first().text()),
                Double.valueOf(cationCoordinateElement.getElementsByTag("y").first().text()),
                Double.valueOf(cationCoordinateElement.getElementsByTag("z").first().text())
        };
        String chainId = element.getElementsByTag("reschain").first().text();
        int residueNumber = Integer.valueOf(element.getElementsByTag("resnr").first().text());
        Optional<Group> optionalCationGroup = structure.select()
                .chainId(chainId)
                .residueNumber(residueNumber)
                .asOptionalGroup();
        if(!optionalCationGroup.isPresent()) {
            return Optional.empty();
        }

        Group cationGroup = optionalCationGroup.get();
        Optional<Atom> optionalCation = cationGroup.atoms()
                .min(Comparator.comparingDouble(atom -> atom.calculate().distanceFast(cationCoordinates)));
        if(!optionalCation.isPresent()) {
            return Optional.empty();
        }
        Atom cation = optionalCation.get();

        List<Atom> piAtoms = element.getElementsByTag("idx")
                .stream()
                .map(Element::text)
                .map(pdbSerial -> selectAtom(structure, pdbSerial))
//                .map(pdbSerial -> selectAtom(filteredAtoms, pdbSerial))
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.toList());
        if(piAtoms.isEmpty()) {
            return Optional.empty();
        }

        Group piGroup = piAtoms.get(0).getParentGroup();
        return Optional.of(new PiCationInteraction(cation,
                piAtoms,
                cationGroup,
                piGroup));
    }

    private void parsePiStackingInteractions(Structure structure, Element interactionElement, List<Interaction> interactions) {
//    private void parsePiStackingInteractions(List<Atom> filteredAtoms, Element interactionElement, List<Interaction> interactions) {
        Elements piStackingInteractions = interactionElement.getElementsByTag("pi_stack");
        piStackingInteractions.stream()
                .map(element -> parsePiStackingInteraction(structure, element))
//                .map(element -> parsePiStackingInteraction(filteredAtoms, element))
                .filter(Optional::isPresent)
                .map(Optional::get)
                .forEach(interactions::add);
    }

    private Optional<PiStackingInteraction> parsePiStackingInteraction(Structure structure, Element element) {
//    private Optional<PiStackingInteraction> parsePiStackingInteraction(List<Atom> filteredAtoms, Element element) {
//        Structure structure = filteredAtoms.get(0).getParentStructure();

        //TODO: potentially half interaction
        Optional<Element> optionalCoordinateElement = element.getElementsByTag("protcoo").stream().findFirst();
        if(!optionalCoordinateElement.isPresent()) {
            return Optional.empty();
        }

        Element coordinateElement = optionalCoordinateElement.get();
        double[] coordinates = new double[] {
                Double.valueOf(coordinateElement.getElementsByTag("x").first().text()),
                Double.valueOf(coordinateElement.getElementsByTag("y").first().text()),
                Double.valueOf(coordinateElement.getElementsByTag("z").first().text())
        };
        String chainId = element.getElementsByTag("reschain").first().text();
        int residueNumber = Integer.valueOf(element.getElementsByTag("resnr").first().text());
        Optional<Group> optionalGroup1 = structure.select()
                .chainId(chainId)
                .residueNumber(residueNumber)
                .asOptionalGroup();
        if(!optionalGroup1.isPresent()) {
            return Optional.empty();
        }

        Group group1 = optionalGroup1.get();
        Optional<Atom> optionalAtom1 = group1.atoms()
                .min(Comparator.comparingDouble(atom -> atom.calculate().distanceFast(coordinates)));
        if(!optionalAtom1.isPresent()) {
            return Optional.empty();
        }
        Atom atom1 = optionalAtom1.get();

        List<Atom> piAtoms = element.getElementsByTag("idx")
                .stream()
                .map(Element::text)
                .map(pdbSerial -> selectAtom(structure, pdbSerial))
//                .map(pdbSerial -> selectAtom(filteredAtoms, pdbSerial))
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.toList());
        if(piAtoms.isEmpty()) {
            return Optional.empty();
        }

        Group piGroup = piAtoms.get(0).getParentGroup();
        return Optional.of(new PiStackingInteraction(atom1,
                piAtoms,
                group1,
                piGroup));
    }

    private void parseSaltBridges(Structure structure, Element interactionElement, List<Interaction> interactions) {
//    private void parseSaltBridges(List<Atom> filteredAtoms, Element interactionElement, List<Interaction> interactions) {
        Elements saltBridges = interactionElement.getElementsByTag("salt_bridge");
        saltBridges.stream()
                .map(element -> parseSaltBridge(structure, element))
//                .map(element -> parseSaltBridge(filteredAtoms, element))
                .filter(Optional::isPresent)
                .map(Optional::get)
                .forEach(interactions::add);
    }

    private Optional<SaltBridge> parseSaltBridge(Structure structure, Element element) {
//    private Optional<SaltBridge> parseSaltBridge(List<Atom> filteredAtoms, Element element) {
//        Structure structure = filteredAtoms.get(0).getParentStructure();

        //TODO: potentially half interaction
        Optional<Element> optionalCoordinateElement = element.getElementsByTag("protcoo").stream().findFirst();
        if(!optionalCoordinateElement.isPresent()) {
            return Optional.empty();
        }

        Element coordinateElement = optionalCoordinateElement.get();
        double[] coordinates = new double[] {
                Double.valueOf(coordinateElement.getElementsByTag("x").first().text()),
                Double.valueOf(coordinateElement.getElementsByTag("y").first().text()),
                Double.valueOf(coordinateElement.getElementsByTag("z").first().text())
        };
        String chainId = element.getElementsByTag("reschain").first().text();
        int residueNumber = Integer.valueOf(element.getElementsByTag("resnr").first().text());
        Optional<Group> optionalGroup1 = structure.select()
                .chainId(chainId)
                .residueNumber(residueNumber)
                .asOptionalGroup();
        if(!optionalGroup1.isPresent()) {
            return Optional.empty();
        }

        Group group1 = optionalGroup1.get();
        Optional<Atom> optionalAtom1 = group1.atoms()
                .min(Comparator.comparingDouble(atom -> atom.calculate().distanceFast(coordinates)));
        if(!optionalAtom1.isPresent()) {
            return Optional.empty();
        }
        Atom atom1 = optionalAtom1.get();

        List<Atom> saltAtoms = element.getElementsByTag("idx")
                .stream()
                .map(Element::text)
                .map(pdbSerial -> selectAtom(structure, pdbSerial))
//                .map(pdbSerial -> selectAtom(filteredAtoms, pdbSerial))
                .filter(Optional::isPresent)
                .map(Optional::get)
                .collect(Collectors.toList());
        if(saltAtoms.isEmpty()) {
            return Optional.empty();
        }

        Group saltGroup = saltAtoms.get(0).getParentGroup();
        return Optional.of(new SaltBridge(atom1,
                saltAtoms,
                group1,
                saltGroup));
    }

    private void parseWaterBridges(Structure structure, Element interactionElement, List<Interaction> interactions) {
//    private void parseWaterBridges(List<Atom> filteredAtoms, Element interactionElement, List<Interaction> interactions) {
        Elements waterBridges = interactionElement.getElementsByTag("water_bridge");
        waterBridges.stream()
                .map(element -> parseWaterBridge(structure, element))
//                .map(element -> parseWaterBridge(filteredAtoms, element))
                .filter(Optional::isPresent)
                .map(Optional::get)
                .forEach(interactions::add);
    }

    private Optional<WaterBridge> parseWaterBridge(Structure structure, Element element) {
        Optional<Atom> optionalDonor = selectAtom(structure, element.getElementsByTag("donor_idx").first().text());
        Optional<Atom> optionalAcceptor = selectAtom(structure, element.getElementsByTag("acceptor_idx").first().text());
//    private Optional<WaterBridge> parseWaterBridge(List<Atom> filteredAtoms, Element element) {
//        Optional<Atom> optionalDonor = selectAtom(filteredAtoms, element.getElementsByTag("donor_idx").first().text());
//        Optional<Atom> optionalAcceptor = selectAtom(filteredAtoms, element.getElementsByTag("acceptor_idx").first().text());
        if(!(optionalDonor.isPresent() && optionalAcceptor.isPresent())) {
            return Optional.empty();
        }
        Atom donor = optionalDonor.get();
        Atom acceptor = optionalAcceptor.get();
        Group donorGroup = donor.getParentGroup();
        Group acceptorGroup = acceptor.getParentGroup();
        return Optional.of(new WaterBridge(donor,
                acceptor,
                donorGroup,
                acceptorGroup));
    }

    private Optional<Atom> selectAtom(Structure structure, String pdbSerialString) {
        int pdbSerial = Integer.valueOf(pdbSerialString);
        return structure.atoms()
                .filter(atom -> atom.getPdbSerial() == pdbSerial)
                .findFirst();
    }

//    private Optional<Atom> selectAtom(List<Atom> filteredAtoms, String pdbSerialString) {
//        return Optional.of(filteredAtoms.get(Integer.valueOf(pdbSerialString) - 1));
//    }

    private String[] composeLigandCommandLineCall(Structure structure,
                                                  Path outputDirectory) throws IOException {
        logger.info("annotating ligand interaction for {}",
                structure.getProteinIdentifier().getPdbId());
        Path inputPath = writeStructureToTemporaryFile(structure);
//       Path inputPath = Paths.get("/home/bittrich/Downloads/3g1h.pdb");
        // base commands - input, output, output format
        return Stream.of("python2",
                getServiceLocation(),
                // results in XML format
                "-x",
                // input file location
                "-f",
                inputPath.toFile().getAbsolutePath(),
                "-o",
                outputDirectory.toFile().getAbsolutePath())
                .toArray(String[]::new);
    }

    private String[] composePolymerCommandLineCall(Structure structure,
                                                   Path outputDirectory,
                                                   String chainId) throws IOException {
        logger.info("annotating polymer interaction for {}",
                structure.getProteinIdentifier().getPdbId());
        Path inputPath = writeStructureToTemporaryFile(structure);
        // base commands - input, output, output format
        return Stream.of("python2",
                getServiceLocation(),
                // results in XML format
                "-x",
                // input file location
                "-f",
                inputPath.toFile().getAbsolutePath(),
                "-o",
                outputDirectory.toFile().getAbsolutePath(),
                // specify polymer mode
                "--intra",
                // specify chain id
                chainId)
                .toArray(String[]::new);
    }
}
