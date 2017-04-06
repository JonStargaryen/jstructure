package de.bioforscher.explorer.membrane.repository;

import de.bioforscher.explorer.model.ExplorerChain;
import org.springframework.data.mongodb.repository.MongoRepository;
import org.springframework.data.mongodb.repository.Query;

import java.io.Serializable;
import java.util.List;

/**
 * Interface to mongoDB for protein chains.
 * Created by bittrich on 2/20/17.
 */
public interface ChainRepository extends MongoRepository<ExplorerChain, Serializable> {
    @Query("{ 'id' : ?0 }")
    List<ExplorerChain> findChainById(String id);
}
