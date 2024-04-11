package gammaspike.operator;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.inference.Operator;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;


@Description("Samples single mutant probability indicator. The hastings ratio is designed so that all subsets of vectors with the same number of 'on' bits are equiprobable")
public class MultinomialSampler extends Operator {
	
	
	final public Input<IntegerParameter> intparameterInput = new Input<>("parameter", "an integer parameter to swap individual values for", Validate.REQUIRED);
	final public Input<RealParameter> pInput = new Input<>("p", "probability of single mutant", Validate.REQUIRED);


	@Override
	public void initAndValidate() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public double proposal() {

		int newValue;
		double qFwd, qBck;
		IntegerParameter parameter = intparameterInput.get();
		double p = pInput.get().getValue();
		
		
		// Sample an index to propose
		final int dim = parameter.getDimension();
        final int pos = Randomizer.nextInt(dim);

        
        // Get current value
        final int oldValue = parameter.getValue(pos);
        if (oldValue < -1 || oldValue > 1) {
        	return Double.NEGATIVE_INFINITY;
        }
        
        /*
        // Resample from prior
		double u = Randomizer.nextFloat();
		
		// Sample old value
		if (oldValue == -1 || oldValue == 1) {
			qBck = Math.log(p/2);
		}else {
			qBck = Math.log(1-p);
		}
		
		
		// Sample new value
		
		if (u < p/2) {
			newValue = -1;
			qFwd = Math.log(p/2);
		}else if (u < p) {
			newValue = 1;
			qFwd = Math.log(p/2);
		}else {
			newValue = 0;
			qFwd = Math.log(1-p);
		}
	

		// Set value
        parameter.setValue(pos, newValue);

		// Return HR
        return qBck - qFwd;
        */
        
        
        
        // Change it to a different value (out of -1, 0, 1)
        if (oldValue == -1) newValue = Randomizer.nextBoolean() ? 0 : 1;
        else if (oldValue == 0) newValue = Randomizer.nextBoolean() ? -1 : 1;
        else newValue = Randomizer.nextBoolean() ? -1 : 0;
        
        
        // How many elements currently have 'value' and 'newValue'
        int nWithValue = 0;
        int nWithNewValue = 0;
        for (int i = 0; i < dim; i++) {
            if (parameter.getValue(i) == oldValue) nWithValue ++;
            if (parameter.getValue(i) == newValue) nWithNewValue ++;
        }
        
        
        // Make the proposal
        parameter.setValue(pos, newValue);
        
        
        // The Hastings ratio is qbck/qfwd = (1/(nWithNewValue+1)) / (1/nWithValue) = nWithValue / (nWithNewValue+1)
        //Log.warning("from " + value + "(" + nWithValue + ") to " + newValue  + "(" + nWithNewValue + ") with HR " + 1.0*nWithValue/(nWithNewValue+1));
        return Math.log(nWithValue) - Math.log(nWithNewValue+1);
        
//
//        double logq = 0.0;
//        if (!value) {
//        	int newValue = Randomizer.nextBoolean() ? -1 : +1;
//            p.setValue(pos, newValue);
//
//            logq = -Math.log((dim - sum) / (sum + 1));
//            //logq += Math.log(2);
//
//        } else {
//            //assert value;
//        	int newValue = 0;
//            p.setValue(pos, newValue);
//            logq = -Math.log(sum / (dim - sum + 1));
//            //logq -= Math.log(2);
//        }
//        return logq;
//	
		
//		double qBck, qFwd;
//		
//		int index = Randomizer.nextInt(parameter.getDimension());
//		
//		

		
	}

}
