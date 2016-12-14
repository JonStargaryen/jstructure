package de.bioforscher.jstructure.alignment.consensus;

import de.bioforscher.jstructure.alignment.AlignmentResult;
import de.bioforscher.jstructure.alignment.svd.SVDSuperimposer;
import de.bioforscher.jstructure.model.BinaryTree;
import de.bioforscher.jstructure.model.Combinatorics;
import de.bioforscher.jstructure.model.Pair;
import de.bioforscher.jstructure.model.TreeNode;
import de.bioforscher.jstructure.model.structure.container.AtomContainer;
import de.bioforscher.jstructure.model.structure.selection.Selection;

import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

/**
 * Merges observations into one tree. This is achieved by searching for the pairs whose distance is minimal and merging
 * (i.e. averaging its coordinates). The process is repeated until no more merging operations can be performed.
 * Created by S on 11.11.2016.
 */
@Deprecated
public class ConsensusTreeComposer extends AbstractConsensusComposer {
    /**
     * all known structures
     */
    private List<AtomContainer> ensembles;
    /**
     * all currently to process combinations of ensembles - a alignment pair stores its rmsd
     */
    private List<AlignmentPair> alignmentPairs;
    private List<TreeNode<AtomContainer>> leaves;
    private List<BinaryTree<AtomContainer>> consensusTrees;
    private AtomContainer consensus;
    private boolean firstEntry;
    private List<? extends AtomContainer> originalContainers;
    private SVDSuperimposer svdSuperimposer;

    public ConsensusTreeComposer() {
        this.alignmentPairs = new ArrayList<>();
        this.consensusTrees = new ArrayList<>();
        this.firstEntry = true;
        this.svdSuperimposer = new SVDSuperimposer(true);
    }

    public void composeConsensusTree(List<? extends AtomContainer> containers) {
        originalContainers = containers;
        ensembles = containers.stream()
            .map(structure -> Selection.on(structure)
                    .cloneElements()
                    .asAtomContainer())
            .collect(Collectors.toList());
        Combinatorics.sequentialPairsOf(ensembles, containers)
                .forEach(pair -> pair.getLeft().setIdentifier(pair.getRight().getIdentifier()));

        logger.debug("processing {} structures", ensembles.size());

        // set up also as tree nodes
        leaves = ensembles.stream()
                .map(TreeNode::new)
                .collect(Collectors.toList());

        // formulate all unordered pairs within the ensemble
        alignmentPairs = Combinatorics.uniquePairsOf(ensembles)
                // compute rmsd and wrap in container
                .map(AlignmentPair::new)
                .collect(Collectors.toList());
        logger.debug("with {} combinations", alignmentPairs.size());

        // iteratively merge pairs
        while(ensembles.size() > 1) {
            logger.debug("current ensemble size: {}", ensembles.size());
            // find and merge closest pair
            AlignmentPair mergedPair = mergeClosestPair();
            consensus = mergedPair.getMergedEntry();
            double rmsd = mergedPair.getAlignmentResult().getRmsd();

            // create tree nodes
            TreeNode<AtomContainer> leftNode;
            TreeNode<AtomContainer> rightNode;
            TreeNode<AtomContainer> consensusNode;

            // both nodes have to be leaves
            if (firstEntry) {
                leftNode = findLeaf(mergedPair.getLeft()).get();
                rightNode = findLeaf(mergedPair.getRight()).get();
                consensusNode = new TreeNode<>(consensus, leftNode, rightNode, rmsd);
                firstEntry = false;
            } else {
                // try to find matching node in existing tree, if node not found in existing tree it has to be a leave
                Optional<TreeNode<AtomContainer>> potentialLeftNode = findNode(mergedPair.getLeft());
                leftNode = potentialLeftNode.orElseGet(() -> findLeaf(mergedPair.getLeft()).get());

                // try to find matching node in existing trees, if node not found in existing tree it has to be a leave
                Optional<TreeNode<AtomContainer>> potentialRightNode = findNode(mergedPair.getRight());
                rightNode = potentialRightNode.orElseGet(() -> findLeaf(mergedPair.getRight()).get());

                consensusNode = new TreeNode<>(consensus, leftNode, rightNode, rmsd);
            }

            // create and set consensus tree
            BinaryTree<AtomContainer> consensusTree = new BinaryTree<>(consensusNode);
            System.out.println("current consensus: " + consensusTree.composeNewickRepresentation());
            consensusTrees.add(consensusTree);
        }
    }

    /**
     * Returns the leave that is equal to the given {@link AtomContainer}.
     * @param observation the {@link AtomContainer} for which a leave should be found
     * @return the corresponding leave
     */
    private Optional<TreeNode<AtomContainer>> findLeaf(final AtomContainer observation) {
        return leaves.stream()
                .filter(leave -> leave.getData().equals(observation))
                .findAny();
    }

    /**
     * Searches in all existing consensus trees to find the node containing the given {@link AtomContainer}.
     * @param observation the {@link AtomContainer} for which a node should be found
     * @return the node containing the given {@link AtomContainer} or null if it was not found
     */
    private Optional<TreeNode<AtomContainer>> findNode(final AtomContainer observation) {
        return consensusTrees.stream()
                .flatMap(BinaryTree::children)
                .filter(node -> node.getData().getIdentifier().equals(observation.getIdentifier()))
                .findAny();
    }

    private AlignmentPair mergeClosestPair() {
        // sort alignmentPairs by rmsd
        final AlignmentPair pairToMerge = alignmentPairs.stream()
                .reduce((a, b) -> a.getAlignmentResult().getRmsd() < b.getAlignmentResult().getRmsd() ? a : b)
                .orElseThrow(() -> new IllegalArgumentException("did not find pair to merge"));

        // remove entries from ensemble
        ensembles.remove(pairToMerge.getLeft());
        ensembles.remove(pairToMerge.getRight());
        // remove entries dealing with now removed entries from alignment pairs
        alignmentPairs.removeIf(alignmentPair -> alignmentPair.describes(pairToMerge));
        logger.debug("pair to merge is {} with rmsd of {}", pairToMerge, pairToMerge.getAlignmentResult().getRmsd());

        // merge entry
        AtomContainer mergedEntry = pairToMerge.getMergedEntry();
        // compute distances for entry and store them
        ensembles.forEach(structure -> alignmentPairs.add(new AlignmentPair(mergedEntry, structure)));
        ensembles.add(mergedEntry);

        return pairToMerge;
    }

    public List<AtomContainer> getAlignedContainers() {
        return originalContainers.stream()
                // align container relative to consensus
                .map(container -> svdSuperimposer.align(consensus, container))
                .map(AlignmentResult::getOriginalQuery)
                .collect(Collectors.toList());
    }

    public BinaryTree<AtomContainer> getConsensusTree() {
        return consensusTrees.get(consensusTrees.size() - 1);
    }

    public AtomContainer getConsensus() {
        return consensus;
    }

    /**
     * Internal helper class to deal with combinations of aligned fragments/containers.
     */
    class AlignmentPair extends Pair<AtomContainer, AtomContainer> {
        private AtomContainer mergedEntry;
        private AlignmentResult alignmentResult;

        AlignmentPair(Pair<? extends AtomContainer, ? extends AtomContainer> alignmentPair) {
            this(alignmentPair.getLeft(), alignmentPair.getRight());
        }

        AlignmentPair(AtomContainer atomContainer1, AtomContainer atomContainer2) {
            super(atomContainer1, atomContainer2);
            this.alignmentResult = svdSuperimposer.align(atomContainer1, atomContainer2);
        }

        AlignmentResult getAlignmentResult() {
            return alignmentResult;
        }

        /**
         * Checks whether this entry shares an element with the alignment pair.
         * @param alignmentPair the entry about to be removed
         * @return true if this entry thus becomes invalid
         */
        boolean describes(AlignmentPair alignmentPair) {
            return contains(alignmentPair.getLeft()) || contains(alignmentPair.getRight());
        }

        AtomContainer getMergedEntry() {
            if(mergedEntry == null) {
                mergedEntry = mergeContainerPair(getLeft(), getRight());
            }

            return mergedEntry;
        }
    }
}
