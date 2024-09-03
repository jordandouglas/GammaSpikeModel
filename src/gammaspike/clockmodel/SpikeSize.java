package gammaspike.clockmodel;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Tree;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;


@Description("multiplies each spike (with mean 1) by the mean spike size")
public class SpikeSize extends CalculationNode implements Function, Loggable {
    final public Input<RealParameter> spikesInput = new Input<>("spikes", "argument to be summed", Validate.REQUIRED);
    final public Input<RealParameter> spikeMeanInput = new Input<>("spikeMean", "mean spike size", Validate.REQUIRED);


    @Override
    public void initAndValidate() {

    }

    @Override
    public int getDimension() {
        return spikesInput.get().getDimension();
    }

    @Override
    public double getArrayValue() {
        return getArrayValue(0);
    }

   

    @Override
    public double getArrayValue(int dim) {
    	return spikesInput.get().getValue(dim) * spikeMeanInput.get().getValue();
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
