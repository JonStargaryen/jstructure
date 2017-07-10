package studies.gmlvq.mutation.preparation;

import de.bioforscher.testutil.TestUtils;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Element;

import java.io.IOException;

/**
 * Fetch a list of all non-redundant PDB chains.
 * Use information from <a href="https://www.ncbi.nlm.nih.gov/Structure/VAST/nrpdbsrv.cgi">https://www.ncbi.nlm.nih.gov/Structure/VAST/nrpdbsrv.cgi</a>,
 * at non-redundancy level: 10e-7.
 * Created by bittrich on 7/10/17.
 */
public class S01_FetchNonRedundantPdbChains {
    public static void main(String[] args) throws IOException {
        Jsoup.parse(TestUtils.getResourceAsInputStream("nrpdb/nrpdb-e70.html"),
                "UTF-8",
                "/")
                .getElementsByTag("TD")
                .stream()
                .map(Element::text)
                .map(string -> string.substring(0, 4).toLowerCase() + "_" + string.substring(4, 5))
                .forEach(System.out::println);
    }
}
