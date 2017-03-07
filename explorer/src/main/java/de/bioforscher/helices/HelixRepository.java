package de.bioforscher.helices;

import de.bioforscher.helices.model.AnnotatedHelix;
import org.springframework.data.mongodb.repository.MongoRepository;
import org.springframework.data.mongodb.repository.Query;

import java.io.Serializable;
import java.util.List;

/**
 * Interface to mongoDB.
 * Created by bittrich on 3/7/17.
 */
public interface HelixRepository extends MongoRepository<AnnotatedHelix, Serializable> {
    @Query("{ 'helixId' : ?0 }")
    List<AnnotatedHelix> findByTheHelixId(String helixId);
}
