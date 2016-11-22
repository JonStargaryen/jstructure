package de.bioforscher.jstructure.model;

import java.util.Optional;
import java.util.stream.Stream;

/**
 * A small binary tree. It defines some hierarchy of {@link TreeNode} objects.
 * Created by S on 15.11.2016.
 */
public class BinaryTree<T> {
    private TreeNode<T> root;

    public BinaryTree(TreeNode<T> root) {
        this.root = root;
    }

    /**
     * Access to the root element of this tree.
     * @return the {@link TreeNode} representing this tree's root element
     */
    public TreeNode<T> getRoot() {
        return root;
    }

    /**
     * Reports the size of this tree.
     * @return this tree's size
     * @see TreeNode#size()
     */
    public int size() {
        return root.size();
    }

    /**
     * Traverses all nodes of this tree and tries to find any node containing the request information.
     * @param data the entry we are looking for
     * @return some node matching this data or {@link Optional#empty()}
     */
    public Optional<TreeNode<T>> findNode(T data) {
        return root.findNode(data);
    }

    /**
     * Stream support for all children elements.
     * @return a select consisting of this element and all children (if any)
     */
    public Stream<TreeNode<T>> children() {
        return root.children();
    }
}
