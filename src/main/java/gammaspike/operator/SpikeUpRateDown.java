package gammaspike.operator;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Scalable;
import beast.base.spec.inference.operator.ScaleOperator;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.util.Randomizer;


@Description("Moves a spike up and the corresponding rate down")
public class SpikeUpRateDown extends ScaleOperator {
	
	final public Input<RealVectorParam<NonNegativeReal>> spikesInput = new Input<>("spikes", "one spike size per branch.", Input.Validate.REQUIRED); 
	

	@Override
	public void initAndValidate() {
		
		Scalable rates = parameterInput.get();
		if (!(rates instanceof RealVectorParam<?>)) {
			throw new IllegalArgumentException("Please ensure the parameter is a RealVectorParam");
		}
		
		
		super.initAndValidate();
		
	}

	@Override
	public double proposal() {

		
		Scalable scalable = parameterInput.get();
		RealVectorParam<?> rates = (RealVectorParam<?>) scalable;
		RealVectorParam<NonNegativeReal> spikes = spikesInput.get();
		
		
		// Sample an index
		final int index = Randomizer.nextInt(rates.size());
		final double scale = getScaler(index, 0);
		
		double r = rates.get(index);
		double s = spikes.get(index);
		
		
		double r_ = r*scale;
		double s_ = s/scale;
		
		
		rates.set(index, r_);
		spikes.set(index, s_);
		
		
		return 0;
	}

}
