package gammaspike.distribution;

import beast.base.evolution.tree.TreeInterface;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;
import junit.framework.TestCase;
import org.junit.Test;

/**
 * @author Robert Haobo Yuan
 */
public class StumpedTreePriorTest extends TestCase {

   @Test
   public void testTreeLikelihoodCalculationOrigin() {
       TreeInterface tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);");

       StumpedTreePrior stumpedTreePrior = new StumpedTreePrior();
       stumpedTreePrior.setInputValue("tree", tree);
       stumpedTreePrior.setInputValue("origin", new RealParameter("10."));
       stumpedTreePrior.setInputValue("conditionOnSampling", false);
       stumpedTreePrior.setInputValue("conditionOnRhoSampling", false);

       stumpedTreePrior.setInputValue("lambda", new RealParameter("2.25"));
       stumpedTreePrior.setInputValue("turnover", new RealParameter("0.46666667"));
       stumpedTreePrior.setInputValue("samplingProportion", new RealParameter("0.3"));
       stumpedTreePrior.setInputValue("rho", new RealParameter("0.0"));

       stumpedTreePrior.setInputValue("ignoreTreePrior", false);
       stumpedTreePrior.setInputValue("ignoreStubPrior", true);

       stumpedTreePrior.initAndValidate();

       // This value was calculated from SA
       assertEquals(-28.7545, stumpedTreePrior.calculateTreeLogLikelihood(tree), 1e-4);

       stumpedTreePrior.setInputValue("tree", tree);
       stumpedTreePrior.setInputValue("origin", new RealParameter("10."));
       stumpedTreePrior.setInputValue("conditionOnSampling", true);
       stumpedTreePrior.setInputValue("conditionOnRhoSampling", false);

       stumpedTreePrior.setInputValue("lambda", new RealParameter("2.25"));
       stumpedTreePrior.setInputValue("turnover", new RealParameter("0.46666667"));
       stumpedTreePrior.setInputValue("samplingProportion", new RealParameter("0.3"));
       stumpedTreePrior.setInputValue("rho", new RealParameter("0.0"));

       stumpedTreePrior.setInputValue("ignoreTreePrior", false);
       stumpedTreePrior.setInputValue("ignoreStubPrior", true);

       stumpedTreePrior.initAndValidate();

       // This value was calculated from SA
       assertEquals(-28.3144, stumpedTreePrior.calculateTreeLogLikelihood(tree), 1e-4);

       stumpedTreePrior.setInputValue("tree", tree);
       stumpedTreePrior.setInputValue("origin", new RealParameter("10."));
       stumpedTreePrior.setInputValue("conditionOnSampling", false);
       stumpedTreePrior.setInputValue("conditionOnRhoSampling", true);

       stumpedTreePrior.setInputValue("lambda", new RealParameter("2.25"));
       stumpedTreePrior.setInputValue("turnover", new RealParameter("0.46666667"));
       stumpedTreePrior.setInputValue("samplingProportion", new RealParameter("0.3"));
       stumpedTreePrior.setInputValue("rho", new RealParameter("0.5"));

       stumpedTreePrior.setInputValue("ignoreTreePrior", false);
       stumpedTreePrior.setInputValue("ignoreStubPrior", true);

       stumpedTreePrior.initAndValidate();

       // This value was calculated from SA
       assertEquals(-30.6927, stumpedTreePrior.calculateTreeLogLikelihood(tree), 1e-4);

   }

    @Test
    public void testTreeLikelihoodCalculationRoot() {
        TreeInterface tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);");

        StumpedTreePrior stumpedTreePrior = new StumpedTreePrior();
        stumpedTreePrior.setInputValue("tree", tree);
        stumpedTreePrior.setInputValue("conditionOnRoot", true);
        stumpedTreePrior.setInputValue("conditionOnSampling", true);
        stumpedTreePrior.setInputValue("conditionOnRhoSampling", false);

        stumpedTreePrior.setInputValue("lambda", new RealParameter("2.25"));
        stumpedTreePrior.setInputValue("turnover", new RealParameter("0.46666667"));
        stumpedTreePrior.setInputValue("samplingProportion", new RealParameter("0.3"));
        stumpedTreePrior.setInputValue("rho", new RealParameter("0.0"));

        stumpedTreePrior.setInputValue("ignoreTreePrior", false);
        stumpedTreePrior.setInputValue("ignoreStubPrior", true);

        stumpedTreePrior.initAndValidate();

        // This value was calculated from SA
        assertEquals(-17.9467, stumpedTreePrior.calculateTreeLogLikelihood(tree), 1e-4);

        stumpedTreePrior.setInputValue("tree", tree);
        stumpedTreePrior.setInputValue("conditionOnRoot", true);
        stumpedTreePrior.setInputValue("conditionOnSampling", false);
        stumpedTreePrior.setInputValue("conditionOnRhoSampling", true);

        stumpedTreePrior.setInputValue("lambda", new RealParameter("2.25"));
        stumpedTreePrior.setInputValue("turnover", new RealParameter("0.46666667"));
        stumpedTreePrior.setInputValue("samplingProportion", new RealParameter("0.3"));
        stumpedTreePrior.setInputValue("rho", new RealParameter("0.5"));

        stumpedTreePrior.setInputValue("ignoreTreePrior", false);
        stumpedTreePrior.setInputValue("ignoreStubPrior", true);

        stumpedTreePrior.initAndValidate();

        // This value was calculated from SA
        assertEquals(-20.1364, stumpedTreePrior.calculateTreeLogLikelihood(tree), 1e-4);

    }

    @Test
    public void testTreeLikelihoodCalculationCanonical() {
        TreeInterface tree = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);");

        StumpedTreePrior stumpedTreePrior = new StumpedTreePrior();
        stumpedTreePrior.setInputValue("tree", tree);
        stumpedTreePrior.setInputValue("origin", new RealParameter("10."));
        stumpedTreePrior.setInputValue("conditionOnSampling", true);
        stumpedTreePrior.setInputValue("conditionOnRhoSampling", false);

        stumpedTreePrior.setInputValue("lambda", new RealParameter("2.25"));
        stumpedTreePrior.setInputValue("mu", new RealParameter("1.05"));
        stumpedTreePrior.setInputValue("psi", new RealParameter("0.45"));
        stumpedTreePrior.setInputValue("rho", new RealParameter("0.0"));

        stumpedTreePrior.setInputValue("ignoreTreePrior", false);
        stumpedTreePrior.setInputValue("ignoreStubPrior", true);

        stumpedTreePrior.initAndValidate();

        // This value was calculated from SA
        assertEquals(-28.3144, stumpedTreePrior.calculateTreeLogLikelihood(tree), 1e-4);

    }

}
