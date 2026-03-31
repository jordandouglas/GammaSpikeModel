package gammaspike.distribution;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Tree;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.inference.distribution.LogNormal;
import beast.base.spec.inference.distribution.TensorDistribution;
import beast.base.spec.type.RealScalar;
import beast.base.spec.type.RealVector;
import beast.base.util.Randomizer;


@Description("Log normal prior but with sampled ancestor leaves removed")
public class BranchRatePrior extends TensorDistribution<RealVector<NonNegativeReal>, Double>  {
	
	final public Input<Tree> treeInput = new Input<>("tree", "tree required for determining if a branch rate is included in model", Input.Validate.REQUIRED); 
	final public Input<RealScalar<PositiveReal>> sigmaInput = new Input<>("sigma", "clock standard deviation (log normal).", Input.Validate.REQUIRED); 
	//final public Input<RealVector<NonNegativeReal>> branchRatesInput = new Input<>("branchRates", "branchRates", Input.Validate.REQUIRED); 
	
	protected LogNormal dist = new LogNormal();
	
	@Override
	public double calculateLogP() {
	
		RealVector<NonNegativeReal> branchRates = paramInput.get();
		logP = calcLogP(branchRates.getElements());
		return logP;
	
	}



	@Override
	public void refresh() {
		// TODO Auto-generated method stub
		
	}



	@Override
	public Double getLowerBoundOfParameter() {
		return 0.0;
	}



	@Override
	public Double getUpperBoundOfParameter() {
		// TODO Auto-generated method stub
		return Double.POSITIVE_INFINITY;
	}



	@Override
	protected double calcLogP(Double... value) {
		return this.calcLogP(Arrays.asList(value));
	}
	
	
    private double calcLogP(List<Double> branchRates) {
    	
    	// Check sigma is positive
        double s = sigmaInput.get().get();
        final double mean = 1;
        if (s <= 0) {
    		return Double.NEGATIVE_INFINITY;
    	}
        double m = Math.log(mean) - (0.5 * s * s);
        
        dist.initByName("M", m, "S", s);
        
        Tree tree = (Tree) treeInput.get();
		int dimension = tree.getNodeCount();
        
		
		double logp = 0;
        for (int nodeNr = 0; nodeNr < dimension; nodeNr ++) {
        	
        	
        	double rate = branchRates.get(nodeNr);
        	if (rate < 0) {
        		logp = Double.NEGATIVE_INFINITY;
        		return logp;
        	}

        	
        	logp += dist.logDensity(rate);
        	
        }
    	
        return logp;
        
    }



	@Override
	public List<Double> sample() {
		
		LogNormal dist = new LogNormal();
		
		List<Double> rates = new ArrayList<Double>();
		
		Tree tree = (Tree) treeInput.get();
		int dimension = tree.getNodeCount();
		
		
		 
        // Check sigma is positive
        double s = sigmaInput.get().get();
        final double mean = 1;
        if (s <= 0) {
        	throw new IllegalArgumentException("Cannot sample branch rates because sigma is non-positive " + s);
    	}
        double m = Math.log(mean) - (0.5 * s * s);
        dist.initByName("M", m, "S", s);
        
        
        for (int nodeNr = 0; nodeNr < dimension; nodeNr ++) {
        	
        	
        	try {
				double spikeOfBranch = dist.inverseCumulativeProbability(Randomizer.nextFloat());
				rates.add(spikeOfBranch);
				
			} catch (Exception e) {
				e.printStackTrace();
				throw new IllegalArgumentException("Unexpected error when sampling from LN(" + m + ", " + s + ")");
			}
        	
        }
        
        
        return rates;
		
		
	}
	


}
