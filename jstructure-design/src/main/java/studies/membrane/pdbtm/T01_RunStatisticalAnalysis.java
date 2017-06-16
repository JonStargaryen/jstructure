package studies.membrane.pdbtm;

import studies.membrane.pdbtm.t01.statistics.*;

/**
 * Runs all tasks of #01.
 * Created by bittrich on 6/9/17.
 */
public class T01_RunStatisticalAnalysis {
    public static void main(String[] args) {
        T011_GeneralStatistics.main(args);
        T012_StatisticsByTopology.main(args);
        T013_StatisticsTransmembraneBySequenceMotifs.main(args);
        T014_StatisticsTransmembraneByInteractionType.main(args);
        T015_StatisticsTransmembraneByInteractingAtoms.main(args);
        T016_StatisticsForInteractingHelices.main(args);
        T017_StatisticsForInteractingHelicesByInteractionType.main(args);
        T018_StatisticsForInteractingHelicesByInteractingAtoms.main(args);
    }
}
