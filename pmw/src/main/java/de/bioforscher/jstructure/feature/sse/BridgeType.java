package de.bioforscher.jstructure.feature.sse;

/**
 * A bridge is formed by two non-overlapping stretches of three residues each
 * (i-1,i,i+1) and (j-1,j,j+1), where i<j.
 * <p>
 * Depending on two basic patterns, a Bridge can be either of type parallel (H
 * bonds in {(i-1,j) and (j,i+1)} OR {(j-1,i) and (i,j-1)}) or antiparallel (H
 * bonds in {(i1,j) and (j,i} OR {(i-1,j+1) and (j-1,i+1)})
 *
 * @author Andreas Prlic
 * @author Aleix Lafita
 *
 */
public enum BridgeType {
    PARALLEL,
    ANTIPARALLEL
}