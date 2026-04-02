package gammaspike.clockmodel;


import java.util.List;

import beast.base.core.Citation;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.util.InputUtil;
import beast.base.spec.domain.NonNegativeInt;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.evolution.branchratemodel.Base;
import beast.base.spec.inference.distribution.ScalarDistribution;
import beast.base.spec.inference.parameter.BoolScalarParam;
import beast.base.spec.inference.parameter.IntVectorParam;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.type.RealScalar;
import beast.base.util.Randomizer;
import gammaspike.tree.Stubs;



@Description("Adds a burst of mutations after speciation events that result in a label change. Each branch has it's own base rate")
@Citation(value =
"Douglas, J., Bouckaert, R., Harris, S.C., Carter Jr, C.W., Wills, P.R. (2025) Evolution is coupled with branching across many granularities of life. Proceedings of the Royal Society Series B 29220250182", DOI = "http://doi.org/10.1098/rspb.2025.0182",
year = 2025, firstAuthorSurname = "Douglas")
public class PunctuatedRelaxedClockModel extends Base implements SpikeModel {
	
	final public Input<Tree> treeInput = new Input<>("tree", "the tree this relaxed clock is associated with.", Input.Validate.REQUIRED);
	final public Input<Stubs> stubsInput = new Input<>("stubs", "stubs of the tree", Input.Validate.OPTIONAL);
	final public Input<IntVectorParam<? extends NonNegativeInt>> nstubsPerBranchInput = new Input<>("nstubsPerBranch", "num stubs per branch.", Input.Validate.OPTIONAL);
	
	final public Input<RealScalarParam<? extends NonNegativeReal>> spikeMeanInput = new Input<>("spikeMean", "mean spike size.", Input.Validate.REQUIRED);
	final public Input<BoolScalarParam> indicatorInput = new Input<>("indicator", "burst size is 0 if this is false", Input.Validate.OPTIONAL);
	final public Input<Boolean> relaxedInput = new Input<>("relaxed", "if false then use strict clock for gradual change", Input.Validate.OPTIONAL);
	final public Input<RealVectorParam<? extends NonNegativeReal>> ratesInput = new Input<>("rates", "the rates associated with nodes in the tree for sampling of individual rates among branches.", Input.Validate.OPTIONAL); 
	final public Input<RealVectorParam<? extends NonNegativeReal>> spikesInput = new Input<>("spikes", "one spike size per branch.", Input.Validate.REQUIRED); 
	final public Input<Boolean> parseFromTreeInput = new Input<>("parseFromTree", "Set to true if initial values are to be loaded from tree metadata.", false); 
	final public Input<ScalarDistribution<RealScalar<? extends NonNegativeReal>, Double>> rateDistInput = new Input<>("distr", "the distribution governing the rates among branches. "
			+ "Must have mean of 1. The clock.rate parameter can be used to change the mean rate.", Input.Validate.OPTIONAL);
	final public Input<Double> initialSpikeSizeInput = new Input<>("initialSpike", "initial value of a spike.", 1.0); 
	final public Input<Boolean> noSpikeOnDatedTipsInput = new Input<>("noSpikeOnDatedTips", "Set to true if dated tips should have a spike of 0.", false); 
	
	
	int nRates;
	double[] ratesArray;
	
	@Override
	public void initAndValidate() {
		
		this.nRates = treeInput.get().getNodeCount();
        this.ratesArray = new double[this.nRates];
        
        if (ratesInput.get() != null && ratesInput.get().size() != this.nRates) {
        	ratesInput.get().setDimension(this.nRates);
        	for (int i = 0; i < this.nRates; i ++) {
        		double val = Randomizer.nextLogNormal(1, 0.5, true);
        		ratesInput.get().set(i, val);
        	}
        }
        
        if (spikesInput.get() != null) {
        	
        	final double initialSpikeSize = initialSpikeSizeInput.get();
        	spikesInput.get().setDimension(this.nRates);
        	for (int i = 0; i < this.nRates; i++) {
        		
        		Node node = treeInput.get().getNode(i);
        		if (node.getParent() != null && node.getParent().isFake()) {
        			spikesInput.get().set(i, 0.0);
        		} else {
        			spikesInput.get().set(i, initialSpikeSize);
        		}
        		
        	}
        }
        
        
        // Initialise rates
        ScalarDistribution<RealScalar<? extends NonNegativeReal>, Double> distribution = rateDistInput.get();
        if (distribution != null) {
        	
	      
	        RealVectorParam<NonNegativeReal> other = new RealVectorParam<NonNegativeReal>();
	        other.setDimension(this.nRates);
	        for (int i = 0; i < other.size(); i ++) {
	        	List<Double> val = distribution.sample();
	        	other.set(i, val.get(0));
	        }
	        
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
        		Object val = node.getMetaData(Stubs.NSTUBS_STR);
        		try {
        			int nstubsOnBranch = (int) ((double)val);
        			nstubsPerBranchInput.get().set(i, nstubsOnBranch);
        			Log.warning("nstubs = " + nstubsOnBranch);
        		} catch(Exception e) {
        			throw new IllegalArgumentException("Cannot parse num stubs from metadata label '" + Stubs.NSTUBS_STR + "' got " + val);
        		}
        		
        		// Parse rate
        		//Log.warning("Parsing rates from tree for node " + i + " using " + ratesInput.get().getID());
        		String var = ratesInput.get().getID();
        		val = node.getMetaData(var);
        		try {
        			double rate = (double) val;
        			ratesInput.get().set(i, rate);
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
        			
        			spikesInput.get().set(i, spike);
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
		
		// If the use-spike indicator is false, burst size is zero
		if (indicatorInput.get() != null && !indicatorInput.get().get()) {
			return 0;
		}
		
		if (noSpikeOnDatedTipsInput.get()) {
			if (node.isLeaf() && node.getHeight() > 0) return 0;
		}
		
		// When estimating the number of stubs directly (as opposed to integrating over),
		// do not count the spike if a) there are no stubs, and b) the sibling of this branch is a sampled ancestor (SA)
		// If stub input is not provided and the sibling is a SA
		if ((stubsInput.get() == null && node.getParent() != null && node.getParent().isFake()) ||
			// If stub input is provided and stub is zero
			// getNStubsOnBranch() used when the number and/or placement of stubs is estimated
			(stubsInput.get() != null && stubsInput.get().getNStubsOnBranch(node.getNr()) == 0 )) {
			
			// If the sibling is a SA, burst size is zero
			if (node.getParent() != null && node.getParent().isFake()) {
				return 0;
			}
			
		}

		double spikeMean = spikeMeanInput.get().get();
		double relativeSpike = spikesInput.get().get(node.getNr());
		return relativeSpike * spikeMean;

	}
	
	public double getBranchRate(Node node) {

		// If rate input not provided, strict clock is used
		if (ratesInput.get() == null) return 1;

		// If the use-relaxed-clock indicator is not provided or true, relaxed clock is used
		if (relaxedInput.get() == null || relaxedInput.get()) {
			return ratesInput.get().get(node.getNr());
		}
		// Otherwise, strict clock is used
		return 1;

	}
	
	@Override
	public double getRateForBranch(Node node) {
		
		// Root has average rate
		double baseRate = meanRateInput.get().get();
		if (node.getLength() <= 0 || node.isDirectAncestor() || node.isRoot()) return baseRate;
		
		double burstSize = getBurstSize(node);
		double relativeBranchRate = getBranchRate(node);
		double branchDistance = node.getLength() * baseRate * relativeBranchRate + burstSize;

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
