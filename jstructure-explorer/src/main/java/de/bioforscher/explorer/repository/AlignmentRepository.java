package de.bioforscher.explorer.repository;

import de.bioforscher.explorer.model.ExplorerAlignment;
import org.springframework.data.mongodb.repository.MongoRepository;
import org.springframework.data.mongodb.repository.Query;

import java.io.Serializable;
import java.util.List;

/**
 * Stores multi-sequence alignments.
 * Created by bittrich on 3/17/17.
 */
public interface AlignmentRepository extends MongoRepository<ExplorerAlignment, Serializable> {
    @Query("{ 'representativeId' : ?0 }")
    List<ExplorerAlignment> findAlignmentByRepresentativeId(String representativeId);
}
