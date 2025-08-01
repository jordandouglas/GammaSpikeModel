package gammaspike.clockmodel;

import org.apache.commons.math.MathException;

import beast.base.core.Citation;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import gammaspike.tree.ForwardTimeSimulatorResub;
import gammaspike.tree.Stubs;



@Description("Adds a burst of mutations after speciation events that result in a label change. Each branch has it's own base rate")
@Citation(value =
"Douglas, J., Bouckaert, R., Harris, S.C., Carter Jr, C.W., Wills, P.R. (2025) Evolution is coupled with branching across many granularities of life. Proceedings of the Royal Society Series B 29220250182", DOI = "http://doi.org/10.1098/rspb.2025.0182",
year = 2025, firstAuthorSurname = "Douglas")
public class PunctuatedRelaxedClockModel extends BranchRateModel.Base implements SpikeModel {
	
	final public Input<Tree> treeInput = new Input<>("tree", "the tree this relaxed clock is associated with.", Input.Validate.REQUIRED);
	final public Input<Stubs> stubsInput = new Input<>("stubs", "stubs of the tree", Input.Validate.OPTIONAL);
	final public Input<IntegerParameter> nstubsPerBranchInput = new Input<>("nstubsPerBranch", "num stubs per brancg.", Input.Validate.OPTIONAL); 
	final public Input<RealParameter> spikeMeanInput = new Input<>("spikeMean", "mean spike size.", Input.Validate.REQUIRED); 
	final public Input<BooleanParameter> indicatorInput = new Input<>("indicator", "burst size is 0 of this is false", Input.Validate.OPTIONAL);
	final public Input<BooleanParameter> relaxedInput = new Input<>("relaxed", "if false then use strict clock", Input.Validate.OPTIONAL);
	final public Input<RealParameter> ratesInput = new Input<>("rates", "the rates associated with nodes in the tree for sampling of individual rates among branches.", Input.Validate.OPTIONAL); 
	final public Input<RealParameter> spikesInput = new Input<>("spikes", "one spike size per branch.", Input.Validate.REQUIRED); 
	final public Input<Boolean> parseFromTreeInput = new Input<>("parseFromTree", "Set to true if initial values are to be loaded from tree metadata.", false); 
	final public Input<ParametricDistribution> rateDistInput = new Input<>("distr", "the distribution governing the rates among branches. "
			+ "Must have mean of 1. The clock.rate parameter can be used to change the mean rate.", Input.Validate.OPTIONAL);
	final public Input<Double> initialSpikeSizeInput = new Input<>("initialSpike", "initial value of a spike.", 1.0); 
	final public Input<Boolean> noSpikeOnDatedTipsInput = new Input<>("noSpikeOnDatedTips", "Set to true if dated tips should have a spike of 0.", false); 
	
	
	int nRates;
	double[] ratesArray;
	
	@Override
	public void initAndValidate() {
		

		
		this.nRates = treeInput.get().getNodeCount();
        this.ratesArray = new double[this.nRates];
        
        
        if (ratesInput.get() != null && ratesInput.get().getDimension() != this.nRates) {
        	ratesInput.get().setDimension(this.nRates);
        	for (int i = 0; i < this.nRates; i ++) {
        		double val = Randomizer.nextLogNormal(1, 0.5, true);
        		ratesInput.get().setValue(i, val);
        	}
        }
        
        if (spikesInput.get() != null) {
        	
        	final double initialSpikeSize = initialSpikeSizeInput.get();
        	spikesInput.get().setDimension(this.nRates);
        	for (int i = 0; i < this.nRates; i++) {
        		
        		Node node = treeInput.get().getNode(i);
        		if (node.getParent() != null && node.getParent().isFake()) {
        			spikesInput.get().setValue(i, 0.0);
        		}else {
        			spikesInput.get().setValue(i, initialSpikeSize);
        		}
        		
        		
        	}
        }
        
        
        
        
        // Initialise rates
        ParametricDistribution distribution = rateDistInput.get();
        if (distribution != null) {
        	
	        Double[][] initialRates0 = null;
			try {
				initialRates0 = distribution.sample(this.nRates);
			} catch (MathException e) {
				e.printStackTrace();
			}
	        Double [] initialRates = new Double[this.nRates];
	        for (int i = 0; i < this.nRates; i++) {
	        	initialRates[i] = initialRates0[i][0];
	        }
	        RealParameter other = new RealParameter(initialRates);
	        ratesInput.get().assignFromWithoutID(other);
        }
        
        
        //Log.warning("nnodes = " + nnodes);
        
        // Parse the initial values from the tree metadata
        if (ratesInput.get() != null && parseFromTreeInput.get()) {
        	
        	Log.warning("Parsing from tree");
        	
        	spikesInput.get().setDimension(this.nRates);
        	ratesInput.get().setDimension(this.nRates);
        	nstubsPerBranchInput.get().setDimension(this.nRates);
        	
        	for (int i = 0; i < this.nRates; i++) {
        		
        		
        		Node node = treeInput.get().getNode(i);
        		
        		// Parse nstubs
        		//Log.warning("Parsing stubs from tree for node " + i + " using " + ForwardTimeSimulatorResub.NSTUBS_STR);
        		Object val = node.getMetaData(ForwardTimeSimulatorResub.NSTUBS_STR);
        		try {
        			int nstubsOnBranch = (int) ((double)val);
        			nstubsPerBranchInput.get().setValue(i, nstubsOnBranch);
        			Log.warning("nstubs = " + nstubsOnBranch);
        		} catch(Exception e) {
        			throw new IllegalArgumentException("Cannot parse num stubs from metadata label '" + ForwardTimeSimulatorResub.NSTUBS_STR + "' got " + val);
        		}
        		
        		// Parse rate
        		//Log.warning("Parsing rates from tree for node " + i + " using " + ratesInput.get().getID());
        		String var = ratesInput.get().getID();
        		val = node.getMetaData(var);
        		try {
        			double rate = (double) val;
        			ratesInput.get().setValue(i, rate);
        			Log.warning("rate = " + rate);
        		} catch(Exception e) {
        			throw new IllegalArgumentException("Cannot parse branch rate from metadata label '" + var + "' got " + val);
        		}
        		
        		
        		// Parse spikes
        		//Log.warning("Parsing spikes from tree for node " + i + " using " + spikesInput.get().getID());
        		var = spikesInput.get().getID();
        		val = node.getMetaData(var);
        		try {
        			double spike = (double) val;
        			
        			
        			// Ensure that all spikes are 0 for sampled ancestors
        			if (node.isDirectAncestor()) spike = 0.0;
        			
        			spikesInput.get().setValue(i, spike);
        			Log.warning("spike = " + spike);
        		} catch(Exception e) {
        			throw new IllegalArgumentException("Cannot parse branch spike from metadata label '" + var + "' got " + val);
        		}
        		
        	}
        	
        	
        	
        }
        


        
	}
	
	
	
	/**
	 * Get the size of a burst (this will be zero if the stubs are estimated directly and the sister is a sampled ancestor)
	 * @param dim
	 * @return
	 */
	public double getBurstSize(int dim) {
		Node node = treeInput.get().getNode(dim);
		return getBurstSize(node);
	}

	
	/**
	 * Get the size of a burst (this will be zero if the stubs are estimated directly and the sister is a sampled ancestor)
	 * @param dim
	 * @return
	 */
	public double getBurstSize(Node node) {
		
		
		if (indicatorInput.get() != null && !indicatorInput.get().getValue()) {
			return 0;
		}
		
		if (noSpikeOnDatedTipsInput.get()) {
			if (node.isLeaf() && node.getHeight() > 0) return 0;
		}
		
		
		// Estimating the number of stubs directly (as opposed to integrating over)
		// Do not count the spike if a) there are no stubs, and b) the sibling of this branch is a sampled ancestor
		if ((stubsInput.get() == null && node.getParent() != null && node.getParent().isFake()) ||
			(stubsInput.get() != null && stubsInput.get().getNStubsOnBranch(node.getNr()) == 0 )) {
			
			// Do not apply the spike if the parent of the branch is a sampled ancestor
			if (node.getParent() != null && node.getParent().isFake()) {
				return 0;
			}
			
		}

		// One spike per branch
		double spikeMean = spikeMeanInput.get().getValue();
		double relativeSpike = spikesInput.get().getValue(node.getNr());
		return relativeSpike * spikeMean;

		
	}
	
	public double getBranchRate(Node node) {
		
		if (ratesInput.get() == null) return 1;
		
		if (relaxedInput.get() == null || relaxedInput.get().getValue()) {
			return ratesInput.get().getValue(node.getNr());
		}
		
		return 1;
		
		
	}
	

	@Override
	public double getRateForBranch(Node node) {
		
		
		
		// Root has average rate
		double baseRate = meanRateInput.get().getArrayValue();
		if (node.getLength() <= 0 || node.isDirectAncestor() || node.isRoot()) return baseRate;
		
		
		
		
		double burstRate = getBurstSize(node);
		double branchRate = getBranchRate(node);
		double branchDistance = node.getLength()*baseRate*branchRate + burstRate;
		
		
		
		// Effective rate takes into account burst and base rate
		double effectiveRate = branchDistance / node.getLength();
		
		//Log.warning(node.getID() + " has burst rate=" + effectiveRate + " b = " + burstRate + " d= " + branchDistance);
		//Log.warning("r=" + effectiveRate);
		

		
		return effectiveRate;
	}




	@Override
    protected boolean requiresRecalculation() {
		

        if (InputUtil.isDirty(spikesInput) || InputUtil.isDirty(meanRateInput) || InputUtil.isDirty(spikeMeanInput) ||
    		InputUtil.isDirty(ratesInput) || InputUtil.isDirty(nstubsPerBranchInput) || InputUtil.isDirty(stubsInput)) {
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
		
		
		for (int i = 0; i < nRates; i ++) {
			Node node = this.treeInput.get().getNode(i);
    		ratesArray[i] = this.getRateForBranch(node);
    		//Log.warning("rate " + i + " = " + ratesArray[i]);
    	}
    	return ratesArray;
	}


	
	

}
