package de.bioforscher.explorer;

import de.bioforscher.explorer.model.ExplorerProtein;
import org.springframework.data.mongodb.repository.MongoRepository;
import org.springframework.data.mongodb.repository.Query;

import java.io.Serializable;
import java.util.List;

/**
 * Interface to mongoDB.
 * Created by bittrich on 2/20/17.
 */
public interface ProteinRepository
         extends MongoRepository<ExplorerProtein, Serializable>
{
    @Query("{ 'name' : ?0 }")
    List<ExplorerProtein> findByTheProteinsName(String name);
}
