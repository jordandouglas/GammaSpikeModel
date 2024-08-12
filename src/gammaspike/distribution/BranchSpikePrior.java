package gammaspike.distribution;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.GammaDistributionImpl;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import gammaspike.tree.Stubs;


@Description("A sum of gamma distributions, one for each spike on a branch")
public class BranchSpikePrior extends Distribution {
	
	final public Input<Stubs> stubsInput = new Input<>("stubs", "stubs of the tree", Input.Validate.OPTIONAL);
	final public Input<IntegerParameter> nstubsInput = new Input<>("nstubs", "num stubs per branch", Input.Validate.XOR, stubsInput);
	final public Input<RealParameter> spikesInput = new Input<>("spikes", "one spike size per branch.", Input.Validate.REQUIRED); 
	final public Input<RealParameter> shapeInput = new Input<>("shape", "shape parameter for the gamma distribution of each spike.", Input.Validate.REQUIRED);
	final public Input<RealParameter> meanInput = new Input<>("mean", "mean (=shape*scale) parameter for the gamma distribution of each spike.", Input.Validate.REQUIRED); 
	
	final public Input<Tree> treeInput = new Input<>("tree", "tree required for setting the spike dimension (if direct sampling)", Input.Validate.OPTIONAL); 
	
	
	org.apache.commons.math.distribution.GammaDistribution gamma = new GammaDistributionImpl(1, 1);
	
	@Override
    public void initAndValidate() {
		 
		 
	}
	
	@Override
	public double calculateLogP() {
        
		logP = 0;
       
        
        // Check shape and scale are positive
        double shape = shapeInput.get().getValue();
        double mean = meanInput.get().getValue();
        if (shape <= 0 || mean <= 0) {
    		logP = Double.NEGATIVE_INFINITY;
    		return logP;
    	}
        double scale = mean / shape;
        
        // Calculate density of the total spike size of each branch, assuming that each node or stub has an iid spike drawn from a Gamma(alpha, beta)
        // This approach integrates across all stub spike sizes, so we don't need to estimate them individually
        Stubs stubs = stubsInput.get();
        for (int nodeNr = 0; nodeNr < spikesInput.get().getDimension(); nodeNr ++) {
        	
        	double spikeOfBranch = spikesInput.get().getValue(nodeNr);
        	if (spikeOfBranch < 0) {
        		logP = Double.NEGATIVE_INFINITY;
        		return logP;
        	}
        	
        	int nstubsOnBranch = stubs == null ? 0 : stubs.getNStubsOnBranch(nodeNr);
        	
        	// Skip sampled ancestors
        	//if (stubs.treeInput.get().getNode(nodeNr).isDirectAncestor()) {
        		//continue;
        	//}
        	
        	//nstubsOnBranch = 0; // TMP
        	
        	double alphaBranch = shape * (nstubsOnBranch + 1); // One spike for the branch, and one per stub
        	gamma = new GammaDistributionImpl(alphaBranch, scale);
        	
        	
        	logP += gamma.logDensity(spikeOfBranch);
        	
        }
        
        // Numerical issue
        if (logP == Double.POSITIVE_INFINITY) {
        	logP = Double.NEGATIVE_INFINITY;
    		return logP;
        }
        
        
        return logP;
    }
	
	@Override
	public List<String> getConditions() {
		List<String> conds = new ArrayList<>();
		conds.add(meanInput.get().getID());
		conds.add(shapeInput.get().getID());
		if (stubsInput.get() != null) conds.add(stubsInput.get().getID());
		if (treeInput.get() != null) conds.add(treeInput.get().getID());
		else conds.add(nstubsInput.get().getID());
		return conds;
	}

	@Override
	public List<String> getArguments() {
		List<String> args = new ArrayList<>();
		args.add(spikesInput.get().getID());
		return args;
	}

	@Override
	public void sample(State state, Random random) {
		
		if (treeInput.get() == null) {
			throw new IllegalArgumentException("Please specify the tree");
		}
		
		if (sampledFlag) return;
		sampledFlag = true;

		// Cause conditional parameters to be sampled
		sampleConditions(state, random);
		
		
		Tree tree = (Tree) treeInput.get();
		int dimension = tree.getNodeCount();
		spikesInput.get().setDimension(dimension);
		//spikesInput.get().setValue(null);
		
	
		
	    // Check shape and scale are positive
        double shape = shapeInput.get().getValue();
        double mean = meanInput.get().getValue();
        if (shape <= 0 || mean <= 0) {
        	throw new IllegalArgumentException("Cannot sample spikes because shape or mean are non-positive " + shape + "  " + mean);
    	}
        double scale = mean / shape;
        
        // Calculate density of the total spike size of each branch, assuming that each node or stub has an iid spike drawn from a Gamma(alpha, beta)
        // This approach integrates across all stub spike sizes, so we don't need to estimate them individually
        Stubs stubs = stubsInput.get();
        for (int nodeNr = 0; nodeNr < spikesInput.get().getDimension(); nodeNr ++) {
        	
        	
        	int nstubsOnBranch = 0;
        	if (nodeNr < nstubsInput.get().getDimension()) {
        		nstubsOnBranch = stubs == null ? nstubsInput.get().getNativeValue(nodeNr) : stubs.getNStubsOnBranch(nodeNr);
        	}
        	
        	
        	
        	//nstubsOnBranch = 0; // TMP
        	
        	
        	double alphaBranch = shape * (nstubsOnBranch + 1); // One spike for the branch, and one per stub
        	gamma = new GammaDistributionImpl(alphaBranch, scale);
        	
        	
        	//Log.warning("Sampling with " + nstubsOnBranch + " stubs");
        	
        	try {
				double spikeOfBranch = gamma.inverseCumulativeProbability(random.nextFloat());
				spikesInput.get().setValue(nodeNr, spikeOfBranch);
				
			} catch (MathException e) {
				e.printStackTrace();
				throw new IllegalArgumentException("Unexpected error when sampling from Gamma(" + shape + ", " + scale + ")");
			}
        	
        }
		
		
		
		
	}
	@Override
    protected boolean requiresRecalculation() {
        return super.requiresRecalculation() || InputUtil.isDirty(stubsInput) || InputUtil.isDirty(spikesInput) | InputUtil.isDirty(shapeInput) | InputUtil.isDirty(meanInput);
    }
	

}
