package de.bioforscher.jstructure.model;

import de.bioforscher.jstructure.model.structure.container.AtomContainer;

import java.util.Optional;
import java.util.stream.Stream;

/**
 * Implementation of a tree node consisting of a combination of a left and a right child tree node.
 * Created by S on 15.11.2016.
 */
public class TreeNode<T> extends Pair<Optional<TreeNode<T>>, Optional<TreeNode<T>>> {
    private final T data;

    public TreeNode(T data) {
        this(data, null, null);
    }

    public TreeNode(T data, TreeNode<T> left, TreeNode<T> right) {
        super(Optional.ofNullable(left), Optional.ofNullable(right));
        this.data = data;
    }

    @SuppressWarnings("unchecked")
    static final TreeNode EMPTY_NODE = new TreeNode(null, null, null) {
        @Override
        public int size() {
            return 0;
        }

        @Override
        public Stream<TreeNode> children() {
            return Stream.empty();
        }
    };

    @SuppressWarnings("unchecked")
    private TreeNode<T> emptyNode() {
        return (TreeNode<T>) EMPTY_NODE;
    }

    /**
     * Reports the size of this node. I.e., the sum of all this children. A node without children has a size of 1.
     * @return the number of elements by traversing down the tree starting from this node
     */
    public int size() {
        // as mentioned the size amounts to 1 and the size of all childrens (if any)
        return 1 + getLeft().orElse(emptyNode()).size() + getRight().orElse(emptyNode()).size();
    }

    /**
     * Access to the data assigned to this tree node. Tree nodes store some hierarchy (its left and right reference) and
     * actual data.
     * @return the stored data
     */
    public T getData() {
        return data;
    }

    /**
     * Stream support for all children elements.
     * @return a select consisting of this element and all children (if any)
     */
    Stream<TreeNode<T>> children() {
        return Stream.concat(Stream.of(this), Stream.concat(getLeft().orElse(emptyNode()).children(),
                getRight().orElse(emptyNode()).children()));
    }

    /**
     * Retrieves a tree node which contains some piece of information.
     * @param data the data the retrieve node should feature
     * @return some node containing the requested data or {@link Optional#empty()}
     */
    public Optional<TreeNode<T>> findNode(final T data) {
        return children()
                .peek(System.out::println)
                .filter(node -> node.data.equals(data))
                .findAny();
    }

    public String composeNewickRepresentation() {
        // return data itself when no children are present
        if(!getLeft().isPresent() || !getRight().isPresent()) {
            return ((AtomContainer) data).getIdentifier();
        }
        return "(" + getLeft().get().composeNewickRepresentation() + "," +
                getRight().get().composeNewickRepresentation() + ")";
    }
}