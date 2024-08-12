package gammaspike.distribution;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;


@Description("Log normal prior but with sampled ancestor leaves removed")
public class BranchRatePrior extends Distribution {
	
	final public Input<Tree> treeInput = new Input<>("tree", "tree required for determining if a branch rate is included in model", Input.Validate.REQUIRED); 
	final public Input<RealParameter> sigmaInput = new Input<>("sigma", "clock standard deviation (log normal).", Input.Validate.REQUIRED); 
	final public Input<RealParameter> branchRatesInput = new Input<>("branchRates", "branchRates", Input.Validate.REQUIRED); 
	
	protected LogNormalImpl dist = new LogNormalImpl(0, 1);
	
	@Override
	public double calculateLogP() {
	
		logP = 0;
       
        // Check sigma is positive
        double s = sigmaInput.get().getValue();
        final double mean = 1;
        if (s <= 0) {
    		logP = Double.NEGATIVE_INFINITY;
    		return logP;
    	}
        double m = Math.log(mean) - (0.5 * s * s);
        dist.setMeanAndStdDev(m, s);
        
        for (int nodeNr = 0; nodeNr < branchRatesInput.get().getDimension(); nodeNr ++) {
        	
        	
        	double rate = branchRatesInput.get().getValue(nodeNr);
        	if (rate < 0) {
        		logP = Double.NEGATIVE_INFINITY;
        		return logP;
        	}
        	
        	// Skip sampled ancestors
        	//if (treeInput.get().getNode(nodeNr).isDirectAncestor()) {
        	//	continue;
        	//}
        	
        	
        	logP += dist.logDensity(rate);
        	
        }
		
        
        return logP;
	
	}
	
	@Override
	public List<String> getArguments() {
		List<String> args = new ArrayList<>();
		args.add(branchRatesInput.get().getID());
		return args;
	}

	@Override
	public List<String> getConditions() {
		List<String> conds = new ArrayList<>();
		conds.add(sigmaInput.get().getID());
		conds.add(treeInput.get().getID());
		return conds;
	}

	@Override
	public void sample(State state, Random random) {
		
		if (sampledFlag) return;
		sampledFlag = true;
		
		
		// Cause conditional parameters to be sampled
		sampleConditions(state, random);
		
		
		Tree tree = (Tree) treeInput.get();
		int dimension = tree.getNodeCount();
		branchRatesInput.get().setDimension(dimension);
		
		
		 
        // Check sigma is positive
        double s = sigmaInput.get().getValue();
        final double mean = 1;
        if (s <= 0) {
        	throw new IllegalArgumentException("Cannot sample branch rates because sigma is non-positive " + s);
    	}
        double m = Math.log(mean) - (0.5 * s * s);
        dist.setMeanAndStdDev(m, s);
        
        
        for (int nodeNr = 0; nodeNr < branchRatesInput.get().getDimension(); nodeNr ++) {
        	
        	
        	try {
				double spikeOfBranch = dist.inverseCumulativeProbability(random.nextFloat());
				branchRatesInput.get().setValue(nodeNr, spikeOfBranch);
				
			} catch (MathException e) {
				e.printStackTrace();
				throw new IllegalArgumentException("Unexpected error when sampling from LN(" + m + ", " + s + ")");
			}
        	
        }
		
        
	}
	
    public class LogNormalImpl implements ContinuousDistribution {
        double m_fMean;
        double m_fStdDev;
        NormalDistributionImpl m_normal = new NormalDistributionImpl(0, 1);

        public LogNormalImpl(double mean, double stdDev) {
            setMeanAndStdDev(mean, stdDev);
        }

        @SuppressWarnings("deprecation")
		void setMeanAndStdDev(double mean, double stdDev) {
            m_fMean = mean;
            m_fStdDev = stdDev;
            m_normal.setMean(mean);
            m_normal.setStandardDeviation(stdDev);
        }

        @Override
        public double cumulativeProbability(double x) throws MathException {
            return m_normal.cumulativeProbability(Math.log(x));
        }

        @Override
        public double cumulativeProbability(double x0, double x1) throws MathException {
            return cumulativeProbability(x1) - cumulativeProbability(x0);
        }

        @Override
        public double inverseCumulativeProbability(double p) throws MathException {
            return Math.exp(m_normal.inverseCumulativeProbability(p));
        }

        @Override
        public double density(double x) {
            if( x <= 0 ) {
                return 0;
            }
            return m_normal.density(Math.log(x)) / x;
        }

        @Override
        public double logDensity(double x) {
            if( x <= 0 ) {
                return  Double.NEGATIVE_INFINITY;
            }
            return m_normal.logDensity(Math.log(x)) - Math.log(x);
        }
    } // class LogNormalImpl

}
