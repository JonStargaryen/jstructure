package design.aggregation;

import java.io.IOException;

import static design.ProteinSource.loadProteins;

/**
 * Parses all data into 1 container.
 *
 * Diary:
 * <pre>
 *     Merely parse proteins and load opm information.
 * </pre>
 * Created by S on 29.10.2016.
 */
public class S03_ExtractTMHelices {
    public static void main(String[] args) throws IOException {
        loadProteins();
    }
}
