package gammaspike.distribution;

import java.io.PrintStream;
import java.util.List;
import java.util.Random;


import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import gammaspike.tree.Stubs;

@Description("Returns posterior density on the stumps in a single lineage, for testing purposes")
public class StumpedTreePriorOneLineage extends Distribution {
	
	final public Input<RealParameter> lambdaInput = new Input<>("lambda", "birth rate lambda", Validate.REQUIRED);
	final public Input<RealParameter> r0Input = new Input<>("r0", "ratio between birth and death rate", Input.Validate.REQUIRED);
	final public Input<RealParameter> tInput = new Input<>("t", "length of branch", Validate.REQUIRED);
	final public Input<Stubs> stubsInput = new Input<>("stubs", "the stubs of this tree", Input.Validate.REQUIRED);

	
	@Override
	public double calculateLogP() {
		
		// Get variables
		double t = tInput.get().getValue();
		double lambda = lambdaInput.get().getValue();
		double mu = lambda / r0Input.get().getValue();
		Stubs stubs = stubsInput.get();
		
		
		logP = calcP(t, lambda, mu, stubs);
		return logP;
	}
	
	
	public double calcP(double t, double lambda, double mu, Stubs stubs) {
		
		double p = 0;
		
		
		
		// Time-dependent process
		int k = stubs.getNStubs();
		double logg = this.getRedRateLogIntegral(lambda, mu, t);
		double g = Math.exp(logg);
		int s = 0;
		for (int i = 0; i < stubs.getStubDimension(); i++) {
			if (!stubs.includeStub(i)) continue;
			s++;
			double stubHeight = stubs.getRelativeTimeOfStub(i) * t;
			p += this.getRedTreeBirthLogRate(stubHeight, lambda, mu);
			//Log.warning("stub " + s + " at height " + stubHeight +  " with rate " + Math.exp(this.getRedTreeBirthLogRate(stubHeight, lambda, mu))); 
		}
		
		
		//p += k*logg;
		p += -g;
		double fact = 1;
		for (int i = 2; i <= k; i ++) {
			fact *= i;
			p += -Math.log(i);
		}
		
		//Log.warning("factorial = " + fact);
		return p;
		
		
	}
	
	

	public double getExpectedNumStubs() {
		
		double t = tInput.get().getValue();
		double lambda = lambdaInput.get().getValue();
		double mu = lambda / r0Input.get().getValue();
		return Math.exp(this.getRedRateLogIntegral(lambda, mu, t));
		
	}
	
	private double getRedRateLogIntegral(double lambda, double mu, double T) {
		double Q = Math.log(lambda-mu) -Math.log(lambda*Math.exp((lambda-mu)*T) - mu);
		double g = 2*Q + 2*lambda*T;
		double logG = Math.log(g);
		return logG;
	}
	
	
	
	private double getRedRateIntegral(double y, double lambda, double mu, double finalTime) {
		double a = Math.log(lambda * Math.exp( (lambda-mu)*(finalTime-y) ) - mu);
		return 2*(a + lambda*y);
	}
	
	
	
	private double getRedTreeBirthLogRate(double timeRemaining, double lambda, double mu) {
		double qt = this.getLogQt(timeRemaining, lambda, mu);
		return Math.log(2) + Math.log(lambda) +  qt;
	}
	
	
	
	/**
	 * Get qt
	 * @param time
	 * @param lambda
	 * @param mu
	 * @return
	 */
	private double getLogQt(double time, double lambda, double mu) {
		double exp = Math.exp((lambda - mu)*time);
		double top = (exp - 1);
		double bottom = lambda*exp - mu;
		double logqt = Math.log(mu) + Math.log(top) - Math.log(bottom);
		return logqt;
	}
	
	
	
    /**
     * Loggable interface implementation follows *
     */
    @Override
    public void init(final PrintStream out) {
        out.print(getID() + "\t" + "stubG\t");
    }

    @Override
    public void log(final long sample, final PrintStream out) {
    	
    	double t = tInput.get().getValue();
		double lambda = lambdaInput.get().getValue();
		double mu = lambda / r0Input.get().getValue();
		double logg = this.getRedRateLogIntegral(lambda, mu, t);
    	double g = Math.exp(logg);
        out.print(getCurrentLogP() + "\t" + g + "\t");
    }

	
	
	
	@Override
    protected boolean requiresRecalculation() {
        return super.requiresRecalculation() || tInput.get().somethingIsDirty()  ||  lambdaInput.get().somethingIsDirty() || r0Input.get().somethingIsDirty() || InputUtil.isDirty(stubsInput);
    }
	
	@Override
	public List<String> getArguments() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<String> getConditions() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void sample(State state, Random random) {
		// TODO Auto-generated method stub
		
	}
	
	
	
	// CDF of this time (not height)
	public double getCDF(double y) {
		
		
		
		double T = tInput.get().getValue();
		double lambda = lambdaInput.get().getValue();
		double mu = lambda / r0Input.get().getValue();
		double g = this.getExpectedNumStubs();
		
		
		if (y <= 0) return 0;
		if (y >= T) return 1;
		
		double p = 1/g;
		p = p * (getRedRateIntegral(y, lambda, mu, T) - getRedRateIntegral(0, lambda, mu, T));
		
		
		if (p < -1e-8 || p > 1 + 1e-8) {
			Log.warning("cdf p = " + p + " for y = " + y + " " + Math.exp(-g) + " " + getRedRateIntegral(y, lambda, mu, T) + " " + getRedRateIntegral(0, lambda, mu, T));
		}
		
		if (p < 0) p = 0;
		if (p > 1) p = 1;
		
		return p;
		
		
	}

	
	// Sample a stub time along [0,T] interval using a numeric approximation of the icdf (bisection method)
	public double sampleTime() {
		
		
		final double EPSILON = 1e-6;
		double T = tInput.get().getValue();
		
		
		double a = 0;
		double b = T;
		double y = (a+b)/2;
		double u = Randomizer.nextDouble();
		double cdf = this.getCDF(y);
		
		// Keep changing y until cdf(y) ~= u
		int ninter = 0;
		while (Math.abs(cdf - u) > EPSILON) {
			
			
			if (cdf < u) {
				
				// Increase y
				a = y;
				
			}else {
				
				// Decrease y
				b = y;
				
			}
			
			y = (a+b)/2;
			
			cdf = this.getCDF(y);
			
			
			ninter++;
			if (ninter > 10000) {
				Log.warning("Warning: cannot solve for u=" + u + ", cdf(y=" + y + ")=" + cdf);
				break;
			}
		}
		
		
		//Log.warning("Warning: successfully solved for u=" + u + ", cdf(y=" + y + ")=" + cdf);
		
		
		// Warning: this is a time not a height!
		return y;
	}


	public double getTotalTime() {
		return tInput.get().getValue();
	}
	

}








