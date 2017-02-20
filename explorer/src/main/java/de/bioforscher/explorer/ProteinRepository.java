package de.bioforscher.explorer;

import de.bioforscher.jstructure.model.structure.Protein;
import org.springframework.data.mongodb.repository.MongoRepository;

import java.io.Serializable;

/**
 * Interface to mongoDB.
 * Created by bittrich on 2/20/17.
 */
public interface ProteinRepository extends MongoRepository<Protein, Serializable> {
}
