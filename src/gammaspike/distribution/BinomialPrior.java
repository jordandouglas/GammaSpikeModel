package gammaspike.distribution;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;

@Description("A prior on a boolean")
public class BinomialPrior extends Distribution {

	
	final public Input<RealParameter> pInput = new Input<>("p", "probabiltiy of seeing a true", Input.Validate.REQUIRED);
	final public Input<BooleanParameter> xInput = new Input<>("x", "", Input.Validate.REQUIRED);
	
	public double calculateLogP() {
        logP = 0;
        
        BooleanParameter x = xInput.get();
        double p = pInput.get().getValue();
        
        if (p < 0 | p > 1) {
        	logP = Double.NEGATIVE_INFINITY;
        	return logP;
        }
        
        int n = x.getDimension();
        int sum = 0;
        for (int i = 0; i < n; i ++) {
        	if (x.getValue(i)) {
        		sum ++;
        	}
        }
        
        // Binomial distribution
        double logFactorialn=0, logFactorial1=0, logFactorial2=0;
        for (int i = 2; i <= n;  i ++) logFactorialn += Math.log(i);
        for (int i = 2; i <= sum; i ++) logFactorial1 += Math.log(i);
        for (int i = 2; i <= (n-sum); i ++) logFactorial2 += Math.log(i);
        
        logP = logFactorialn - (logFactorial1 + logFactorial2) + sum*Math.log(p) + (n-sum)*Math.log(1-p);
        
        return logP;
	}
	
	@Override
	public List<String> getArguments() {
		List<String> conds = new ArrayList<>();
		conds.add(xInput.get().getID());
		return conds;
	}

	@Override
	public List<String> getConditions() {
		List<String> conds = new ArrayList<>();
		conds.add(pInput.get().getID());
		return conds;
	}

	@Override
	public void sample(State state, Random random) {
		
		
		if (this.sampledFlag) return;
		this.sampleConditions(state, random);
		BooleanParameter x = xInput.get();
        double p = pInput.get().getValue();
        int n = x.getDimension();
        
        
		for (int i = 0; i < n; i ++) {
			if (random.nextFloat() < p) {
				x.setValue(i, true);
			}else {
				x.setValue(i, false);
			}
		}
		
		
		this.sampledFlag = true;
	}

}
