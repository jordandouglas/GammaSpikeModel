package gammaspike.operator;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.operator.kernel.BactrianScaleOperator;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;


@Description("Moves a spike up and the corresponding rate down")
public class SpikeUpRateDown extends BactrianScaleOperator {
	
	final public Input<RealParameter> spikesInput = new Input<>("spikes", "one spike size per branch.", Input.Validate.REQUIRED); 
	

	@Override
	public void initAndValidate() {
		super.initAndValidate();
		
	}

	@Override
	public double proposal() {

		
		RealParameter rates = parameterInput.get();
		RealParameter spikes = spikesInput.get();
		
		
		// Sample an index
		final int index = Randomizer.nextInt(rates.getDimension());
		final double scale = getScaler(index, 0);
		
		double r = rates.getValue(index);
		double s = spikes.getValue(index);
		
		
		double r_ = r*scale;
		double s_ = s/scale;
		
		
		rates.setValue(index, r_);
		spikes.setValue(index, s_);
		
		
		return 0;
	}

}
