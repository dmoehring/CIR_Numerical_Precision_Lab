package lab.cir;

import lab.cir.config.ExperimentConfig;
import lab.cir.model.BermudanCallOption;
import lab.cir.model.CashflowSchedule;
import lab.cir.model.CirParameters;

import java.util.Properties;

/** Shared test fixtures: the standard case of the specification (section 19.1). */
public final class TestFixtures {

    private TestFixtures() {
    }

    public static CirParameters model() {
        return CirParameters.of("0.50", "0.04", "0.10", "0.03");
    }

    public static CashflowSchedule bond() {
        return CashflowSchedule.couponBond("100.00", "0.04", "5.0", 1);
    }

    public static BermudanCallOption bermudan() {
        return BermudanCallOption.of(bond(), "100.00", "1.0", "2.0", "3.0", "4.0");
    }

    public static Properties properties(int paths) {
        Properties p = new Properties();
        p.setProperty("cir.initialRate", "0.03");
        p.setProperty("cir.meanReversionSpeed", "0.50");
        p.setProperty("cir.longTermMean", "0.04");
        p.setProperty("cir.volatility", "0.10");
        p.setProperty("bond.notional", "100.00");
        p.setProperty("bond.couponRate", "0.04");
        p.setProperty("bond.maturity", "5.0");
        p.setProperty("bond.couponFrequency", "1");
        p.setProperty("option.strike", "100.00");
        p.setProperty("option.exerciseDates", "1.0,2.0,3.0,4.0");
        p.setProperty("simulation.stepsPerYear", "12");
        p.setProperty("simulation.paths", String.valueOf(paths));
        p.setProperty("simulation.seed", "42");
        p.setProperty("valuation.paths", String.valueOf(paths));
        p.setProperty("valuation.seed", "4242");
        p.setProperty("cluster.count", "5");
        p.setProperty("cluster.minPathsPerCluster", "20");
        p.setProperty("decimal.precision", "34");
        p.setProperty("output.csv", "false");
        return p;
    }

    public static ExperimentConfig config(int paths) {
        return ExperimentConfig.of(properties(paths));
    }
}
