package de.bioforscher.jstructure.membrane.division;

public class A02_MoveDynaMineResults {
//    public static void main(String[] args) {
//        MembraneConstants.list(Paths.get("/home/bittrich/Downloads/dynamine/dynamine/"))
//                .filter(Files::isDirectory)
//                .filter(path -> path.toFile().getName().contains("_"))
//                .map(path -> MembraneConstants.list(path).findFirst())
//                .filter(Optional::isPresent)
//                .map(Optional::get)
//                .map(path -> MembraneConstants.list(path).filter(p -> p.toFile().getName().endsWith(".pred")).findFirst())
//                .filter(Optional::isPresent)
//                .map(Optional::get)
//                .forEach(path -> MembraneConstants.move(path, MembraneConstants.FOLDING_CORES_DIRECTORY
//                            .resolve("division")
//                            .resolve("dynamine")
//                            .resolve(path.toFile().getName())));
//    }

//    public static void main(String[] args) {
//        MembraneConstants.list(MembraneConstants.DIVISION_DIRECTORY.resolve("dynamine"))
//                .forEach(path -> {
//                    String[] split = path.toFile().getName().split("_");
//                    String pdbId = split[0];
//
//                    if(Files.exists(MembraneConstants.PDBTM_NR_ALPHA_DATASET_PDB_DIRECTORY.resolve(pdbId + ".pdb"))) {
//                        MembraneConstants.move(path, MembraneConstants.PDBTM_NR_ALPHA_DATASET_DIRECTORY.resolve("dynamine").resolve(path.toFile().getName()));
//                    }
//                    if(Files.exists(MembraneConstants.PDBTM_NR_BETA_DATASET_PDB_DIRECTORY.resolve(pdbId + ".pdb"))) {
//                        MembraneConstants.move(path, MembraneConstants.PDBTM_NR_BETA_DATASET_DIRECTORY.resolve("dynamine").resolve(path.toFile().getName()));
//                    }
//                });
//    }
}
