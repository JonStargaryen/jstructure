package de.bioforscher.explorer.membrane.repository;

import de.bioforscher.explorer.membrane.model.MultiSequenceAlignment;
import org.springframework.data.mongodb.repository.MongoRepository;
import org.springframework.data.mongodb.repository.Query;

import java.io.Serializable;
import java.util.List;

/**
 * Stores multi-sequence alignments.
 * Created by bittrich on 3/17/17.
 */
public interface AlignmentRepository extends MongoRepository<MultiSequenceAlignment, Serializable> {
    @Query("{ 'representativeChainId' : ?0 }")
    List<MultiSequenceAlignment> findAlignmentByRepresentativeChainId(String representativeChainId);
}
