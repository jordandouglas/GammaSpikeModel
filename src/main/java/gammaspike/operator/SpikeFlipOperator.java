package gammaspike.operator;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.inference.StateNode;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.util.Randomizer;


@Description("Flips a spike from zero to non zero, or vice versa")
public class SpikeFlipOperator extends Operator {
	
	final public Input<RealVectorParam<NonNegativeReal>> spikesInput = new Input<>("spikes", "one spike size per branch.", Input.Validate.REQUIRED); 
	final public Input<Double> spikeMeanInput = new Input<>("spikeMean", "mean of the exponential distribuution for making spikes.", 0.005); 
	

	@Override
	public void initAndValidate() {
		
		
	}

	@Override
	public double proposal() {

		
		double spikeMean = spikeMeanInput.get();
		double lambda = 1.0/spikeMean;
		RealVectorParam<NonNegativeReal> spikes = spikesInput.get();
		
		
		// Sample an index
		final int index = Randomizer.nextInt(spikes.size());
		double sOld = spikes.get(index);
		double sNew = 0;
		double pFwd, pRev;
		
		
		// Add a spike
		if (sOld == 0) {
			
			
			sNew = Randomizer.nextExponential(lambda);
			pRev = 0;
			pFwd = Math.log(lambda) - lambda*sNew; // Exponential density
			
		}
		
		// Delete the spike
		else {
			
			sNew = 0;
			pFwd = 0;
			pRev = Math.log(lambda) - lambda*sOld; // Exponential density
			
		}
		
		spikes.set(index, sNew);
		
		return pRev - pFwd;
	}
	
	
	@Override
	public List<StateNode> listStateNodes() {
		final List<StateNode> list = new ArrayList<>();
		list.add(spikesInput.get());
		return list;
   }

	
}
