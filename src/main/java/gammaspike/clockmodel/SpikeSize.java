package gammaspike.clockmodel;

import java.io.PrintStream;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.Input.Validate;
import beast.base.inference.CalculationNode;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;


@Description("multiplies each spike (with mean 1) by the mean spike size")
public class SpikeSize extends CalculationNode implements Function, Loggable {
    final public Input<RealVectorParam<NonNegativeReal>> spikesInput = new Input<>("spikes", "argument to be summed", Validate.REQUIRED);
    final public Input<RealScalarParam<NonNegativeReal>> spikeMeanInput = new Input<>("spikeMean", "mean spike size", Validate.REQUIRED);
    
    final public Input<PunctuatedRelaxedClockModel> clockModelInput = 
    		new Input<>("clockModel", "clockModel is only required the spikes are estimated directly (not the default behaviour) and when sampled ancestors are used", Validate.OPTIONAL);
    


    @Override
    public void initAndValidate() {

    }

    @Override
    public int getDimension() {
        return spikesInput.get().size();
    }

    @Override
    public double getArrayValue() {
        return getArrayValue(0);
    }

   

    @Override
    public double getArrayValue(int dim) {
    	
    	if (clockModelInput.get() != null) {
    		return clockModelInput.get().getBurstSize(dim);
    	}else {
    		return spikesInput.get().get(dim) * spikeMeanInput.get().get();
    	}
    	
    	
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
    		out.print(this.getArrayValue(i) + "\t");
    	}
    }

    @Override
    public void close(PrintStream out) {
        // nothing to do
    }

} // class Sum
