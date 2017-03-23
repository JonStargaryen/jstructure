package de.bioforscher.explorer.membrane.repository;

import de.bioforscher.explorer.membrane.model.ExplorerAlignment;
import org.springframework.data.mongodb.repository.MongoRepository;
import org.springframework.data.mongodb.repository.Query;

import java.io.Serializable;
import java.util.List;

/**
 * Stores multi-sequence alignments.
 * Created by bittrich on 3/17/17.
 */
public interface AlignmentRepository extends MongoRepository<ExplorerAlignment, Serializable> {
    @Query("{ 'representativeChainId' : ?0 }")
    List<ExplorerAlignment> findAlignmentByRepresentativeChainId(String representativeChainId);
}
