package gammaspike.clockmodel;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.Input.Validate;
import beast.base.inference.CalculationNode;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.type.RealVector;


@Description("multiplies each spike (with mean 1) by the mean spike size")
public class SpikeSize extends CalculationNode implements RealVector<NonNegativeReal>, Loggable {
    final public Input<RealVectorParam<NonNegativeReal>> spikesInput = new Input<>("spikes", "argument to be summed", Validate.REQUIRED);
    final public Input<RealScalarParam<NonNegativeReal>> spikeMeanInput = new Input<>("spikeMean", "mean spike size", Validate.REQUIRED);
    
    final public Input<PunctuatedRelaxedClockModel> clockModelInput = 
    		new Input<>("clockModel", "clockModel is only required the spikes are estimated directly (not the default behaviour) and when sampled ancestors are used", Validate.OPTIONAL);
    


    @Override
    public void initAndValidate() {

    }

    public int getDimension() {
        return spikesInput.get().size();
    }

   

    /**
     * Loggable interface implementation follows
     */
    @Override
    public void init(PrintStream out) {
    	for (int i = 0; i < this.getDimension(); i ++) {
    		String id = this.getID();
    		if (id == null || id.equals("")) id = "weightedSpike";
    		out.print(id + "." + i + "\t");
    	}
    }

    @Override
    public void log(long sampleNr, PrintStream out) {
    	for (int i = 0; i < this.getDimension(); i ++) {
    		out.print(this.get(i) + "\t");
    	}
    }

    @Override
    public void close(PrintStream out) {
        // nothing to do
    }

	@Override
	public List<Double> getElements() {
		List<Double> elements = new ArrayList<>();
		for (int i = 0; i < spikesInput.get().size(); i ++) {
			elements.add(get(i));
		}
		return elements;
	}

	@Override
	public NonNegativeReal getDomain() {
		return NonNegativeReal.INSTANCE;
	}

	@Override
	public double get(int i) {
		if (clockModelInput.get() != null) {
    		return clockModelInput.get().getBurstSize(i);
    	}else {
    		return spikesInput.get().get(i) * spikeMeanInput.get().get();
    	}
	}

} 
