package de.bioforscher.explorer.model;

import java.util.List;

/**
 * The data structure of a mutated position.
 * Created by bittrich on 4/20/17.
 */
public class ExplorerMutation {
    private String id, pdb;
    private List<ExplorerGroup> groups;
    private List<ExplorerLigand> ligands;
}
