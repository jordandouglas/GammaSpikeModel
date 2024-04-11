package gammaspike.clockmodel;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import gammaspike.flabel.Flabel;
import gammaspike.tree.ForwardTimeSimulatorResub;
import gammaspike.tree.Stubs;



@Description("Adds a burst of mutations after speciation events that result in a label change. Each branch has it's own base rate")
public class PunctuatedRelaxedClockModel extends BranchRateModel.Base implements SpikeModel {
	
	final public Input<Tree> treeInput = new Input<>("tree", "the tree this relaxed clock is associated with.", Input.Validate.REQUIRED);
	
	final public Input<Flabel> flabelsInput = new Input<>("flabel", "leaf and internal labels", Input.Validate.OPTIONAL);
	final public Input<Stubs> stubsInput = new Input<>("stubs", "stubs of the tree", Input.Validate.OPTIONAL);
	final public Input<IntegerParameter> nstubsPerBranchInput = new Input<>("nstubsPerBranch", "num stubs per brancg.", Input.Validate.OPTIONAL); 
	
	
	
	final public Input<BooleanParameter> indicatorInput = new Input<>("indicator", "burst size is 0 of this is false", Input.Validate.OPTIONAL);
	final public Input<BooleanParameter> relaxedInput = new Input<>("relaxed", "if false then use strict clock", Input.Validate.OPTIONAL);
	
	final public Input<RealParameter> ratesInput = new Input<>("rates", "the rates associated with nodes in the tree for sampling of individual rates among branches.", Input.Validate.REQUIRED); 
	
	final public Input<RealParameter> burstSizeInput = new Input<>("burstSize", "the additive clock rate after a label change.", Input.Validate.OPTIONAL); 
	final public Input<RealParameter> spikesInput = new Input<>("spikes", "one spike size per branch.", Input.Validate.XOR, burstSizeInput); 
	
	final public Input<Boolean> parseFromTreeInput = new Input<>("parseFromTree", "Set to true if initial values are to be loaded from tree metadata.", false); 
	
	
	
	
	int nnodes;
	double[] ratesArray;
	
	@Override
	public void initAndValidate() {
		
		if (flabelsInput.get() == null && stubsInput.get() == null && nstubsPerBranchInput.get() == null) {
			throw new IllegalArgumentException("Please specify one of: flabels, stubs, nstubsPerBranch");
		}
		
		
		this.nnodes = treeInput.get().getNodeCount();
        this.ratesArray = new double[this.nnodes];
        
        if (ratesInput.get().getDimension() != this.nnodes) {
        	ratesInput.get().setDimension(this.nnodes);
        	for (int i = 0; i < this.nnodes; i ++) {
        		//rateInput.get().setValue(i, 1.0);
        	}
        }
        
        if (spikesInput.get() != null) {
        	
        	final double initialSpikeSize = 0.1 * treeInput.get().getRoot().getHeight();
        	spikesInput.get().setDimension(this.nnodes);
        	for (int i = 0; i < this.nnodes; i++) {
        		spikesInput.get().setValue(i, initialSpikeSize);
        	}
        }
        
        
        //Log.warning("nnodes = " + nnodes);
        
        // Parse the initial values from the tree metadata
        if (parseFromTreeInput.get()) {
        	
        	spikesInput.get().setDimension(this.nnodes);
        	ratesInput.get().setDimension(this.nnodes);
        	nstubsPerBranchInput.get().setDimension(this.nnodes);
        	for (int i = 0; i < this.nnodes; i++) {
        		
        		
        		Node node = treeInput.get().getNode(i);
        		
        		// Parse nstubs
        		Object val = node.getMetaData(ForwardTimeSimulatorResub.NSTUBS_STR);
        		try {
        			int nstubsOnBranch = (int) ((double)val);
        			nstubsPerBranchInput.get().setValue(i, nstubsOnBranch);
        			Log.warning("nstubs = " + nstubsOnBranch);
        		} catch(Exception e) {
        			throw new IllegalArgumentException("Cannot parse num stubs from metadata label '" + ForwardTimeSimulatorResub.NSTUBS_STR + "' got " + val);
        		}
        		
        		// Parse rate
        		String var = ratesInput.get().getID();
        		val = node.getMetaData(var);
        		try {
        			double rate = (double) val;
        			ratesInput.get().setValue(i, rate);
        			Log.warning("rate = " + rate);
        		} catch(Exception e) {
        			throw new IllegalArgumentException("Cannot parse branch rate from metadata label '" + var + "' got " + val);
        		}
        		
        		
        		// Parse rate
        		var = spikesInput.get().getID();
        		val = node.getMetaData(var);
        		try {
        			double spike = (double) val;
        			spikesInput.get().setValue(i, spike);
        			Log.warning("spike = " + spike);
        		} catch(Exception e) {
        			throw new IllegalArgumentException("Cannot parse branch spike from metadata label '" + var + "' got " + val);
        		}
        		
        	}
        	
        	
        	
        }
        
        	
       
        	
        

        
	}

	
	public double getBurstSize(Node node, int nbursts) {
		if (indicatorInput.get() != null && !indicatorInput.get().getValue()) {
			return 0;
		}
		
		// One spike for the whole tree
		if (burstSizeInput.get() != null) {
			return nbursts * burstSizeInput.get().getValue();
		}
		
		// One spike per branch
		else {
			return spikesInput.get().getValue(node.getNr());
		}

		
	}
	
	public double getBranchRate(Node node) {
		
		if (relaxedInput.get() == null || relaxedInput.get().getValue()) {
			return ratesInput.get().getValue(node.getNr());
		}
		
		return 1;
		
		
	}
	

	@Override
	public double getRateForBranch(Node node) {
		
		
		// Root has average rate
		double baseRate = meanRateInput.get().getArrayValue();
		if (node.isRoot()) return baseRate;
		
		// Was there a change in function along this branch?
		
		//flabels.update();
		
		int numberOfSBursts = getNumStubsOnBranch(node); 
		//int numberOfDBursts = flabels.getNumberOfDoubleBursts(node); // tmp: Ignore double bursts
		
		
		double burstRate = getBurstSize(node, numberOfSBursts);
				
		double branchRate = getBranchRate(node);
		double branchDistance = node.getLength()*baseRate*branchRate + burstRate;
		
		
		
		// Effective rate takes into account burst and base rate
		double effectiveRate = branchDistance / node.getLength();
		
		//Log.warning(node.getID() + " has burst rate=" + effectiveRate + " b = " + burstRate + " d= " + branchDistance);
		
		return effectiveRate;
	}

	
	private int getNumStubsOnBranch(Node node) {
		
		Flabel flabels = flabelsInput.get();
		if (flabels != null) {
			return flabels.getNumberOfSingleBursts(node); 
		}
		
		Stubs stubs = stubsInput.get();
		if (stubs != null) {
			return stubs.getNStubsOnBranch(node.getNr());
		}
		
		return nstubsPerBranchInput.get().getNativeValue(node.getNr());
	}


	@Override
    protected boolean requiresRecalculation() {

        if (InputUtil.isDirty(burstSizeInput) || InputUtil.isDirty(spikesInput) || InputUtil.isDirty(meanRateInput) || 
    		InputUtil.isDirty(ratesInput) || InputUtil.isDirty(nstubsPerBranchInput) || InputUtil.isDirty(flabelsInput) || InputUtil.isDirty(stubsInput)) {
       	 	return true;
        }
        
        if (relaxedInput.get() != null && InputUtil.isDirty(relaxedInput)) {
        	return true;
        }
        
        if (indicatorInput.get() != null && InputUtil.isDirty(indicatorInput)) {
        	return true;
        }

        return false;
    }


	@Override
	public double[] getRatesArray() {
		for (int i = 0; i < nnodes; i ++) {
			Node node = this.treeInput.get().getNode(i);
    		ratesArray[i] = this.getRateForBranch(node);
    		//Log.warning("rate " + i + " = " + ratesArray[i]);
    	}
    	return ratesArray;
	}
	

}
