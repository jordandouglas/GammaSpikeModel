package gammaspike.distribution;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.spec.domain.UnitInterval;
import beast.base.spec.inference.distribution.TensorDistribution;
import beast.base.spec.type.BoolVector;
import beast.base.spec.type.RealScalar;
import beast.base.util.Randomizer;

@Description("A prior on a boolean")
public class BinomialPrior extends TensorDistribution<BoolVector, Boolean> {

	
	final public Input<RealScalar<UnitInterval>> pInput = new Input<>("p", "probabiltiy of seeing a true", Input.Validate.REQUIRED);
	
	
	@Override
	public double calculateLogP() {
		BoolVector x = paramInput.get();
		logP = calcLogP(x.getElements());
		return logP;
	}




	@Override
	public void refresh() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public Boolean getLowerBoundOfParameter() {
		throw new IllegalStateException(getClass().getName() + " does not support lower bounds.");
	}

	@Override
	public Boolean getUpperBoundOfParameter() {
		 throw new IllegalStateException(getClass().getName() + " does not support upper bounds.");
	}

	@Override
	protected double calcLogP(Boolean... value) {
		return this.calcLogP(Arrays.asList(value));
    }

    private double calcLogP(List<Boolean> x) {
    	
        double p = pInput.get().get();
        
        if (p < 0 | p > 1) {
        	return Double.NEGATIVE_INFINITY;
        }
        
        int n = x.size();
        int sum = 0;
        for (int i = 0; i < n; i ++) {
        	if (x.get(i)) {
        		sum ++;
        	}
        }
        
        // Binomial distribution
        double logFactorialn=0, logFactorial1=0, logFactorial2=0;
        for (int i = 2; i <= n;  i ++) logFactorialn += Math.log(i);
        for (int i = 2; i <= sum; i ++) logFactorial1 += Math.log(i);
        for (int i = 2; i <= (n-sum); i ++) logFactorial2 += Math.log(i);
        
        double logp = logFactorialn - (logFactorial1 + logFactorial2) + sum*Math.log(p) + (n-sum)*Math.log(1-p);
        return logp;
        
    	
    }
	
	

	@Override
	public List<Boolean> sample() {
		
		List<Boolean> vals = new ArrayList<Boolean>();
		
        double p = pInput.get().get();
        int n = paramInput.get().size();
        
		for (int i = 0; i < n; i ++) {
			if (Randomizer.nextFloat() < p) {
				vals.add(true);
			}else {
				vals.add(false);
			}
		}
		
		return vals;
	}


}
