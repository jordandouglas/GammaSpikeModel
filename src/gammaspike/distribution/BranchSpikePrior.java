package gammaspike.distribution;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;


import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.GammaDistributionImpl;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import gammaspike.tree.Stubs;


@Description("A sum of gamma distributions, one for each spike on a branch")
public class BranchSpikePrior extends Distribution {
	
	final public Input<Stubs> stubsInput = new Input<>("stubs", "stubs of the tree", Input.Validate.OPTIONAL);
	final public Input<IntegerParameter> nstubsInput = new Input<>("nstubs", "num stubs per branch", Input.Validate.OPTIONAL);
	final public Input<RealParameter> spikesInput = new Input<>("spikes", "one spike size per branch.", Input.Validate.REQUIRED); 
	final public Input<RealParameter> shapeInput = new Input<>("shape", "shape parameter for the gamma distribution of each spike.", Input.Validate.REQUIRED);
	final public Input<RealParameter> meanInput = new Input<>("mean", "mean (=shape*scale) parameter for the gamma distribution of each spike.", Input.Validate.OPTIONAL); 
	
	final public Input<Tree> treeInput = new Input<>("tree", "tree required for setting the spike dimension (if direct sampling)", Input.Validate.OPTIONAL); 
	final public Input<BooleanParameter> indicatorInput = new Input<>("indicator", "burst size is 0 of this is false", Input.Validate.OPTIONAL);
	
	//final public Input<Boolean> meanGamma1Input = new Input<>("meanGamma1", "if true, then spikes are divided by 'mean' to avoid smaller values of 'mean' from having a larger density", true);
	
	
	// If there are too many stubs on a branch (eg. during mixing) then the gamma distribution shape is large, which causes
	// instbilities
	final double MAX_CUM_SUM = 0.999;
	
	org.apache.commons.math.distribution.GammaDistribution gamma = new GammaDistributionImpl(1, 1);
	
	@Override
    public void initAndValidate() {
		 
		if (stubsInput.get() != null) {
			stubsInput.get().setBranchSpikePrior(this);
		}
		 
	}
	
	@Override
	public double calculateLogP() {
        
		logP = 0;
       
        
        // Check shape and scale are positive
        double shape = shapeInput.get().getValue();
        double mean = 1; //meanInput.get().getValue();
        if (shape <= 0 || mean <= 0) {
    		logP = Double.NEGATIVE_INFINITY;
    		return logP;
    	}
        double scale = mean / shape;
        
        //if (meanGamma1Input.get()) {
    		//scale = 1/shape;
    	//}
        
        
//        double xxx = 1e-30;
//        for (int i = 0; i < 30; i ++) {
//        	xxx *= 10;
//        	gamma = new GammaDistributionImpl(2, 0.5);
//        	double p = gamma.logDensity(xxx);
//        	Log.warning(xxx + " " + p);
//        	
//        }
        
        // Calculate density of the total spike size of each branch, assuming that each node or stub has an iid spike drawn from a Gamma(alpha, beta)
        // This approach integrates across all stub spike sizes, so we don't need to estimate them individually
        Stubs stubs = stubsInput.get();
        for (int nodeNr = 0; nodeNr < spikesInput.get().getDimension(); nodeNr ++) {
        	
        	double spikeOfBranch = spikesInput.get().getValue(nodeNr);
        	
        	
        	
        	if (spikeOfBranch < 0) {
        		logP = Double.NEGATIVE_INFINITY;
        		return logP;
        	}
        	

        	
        	if (stubs == null || stubs.estimateStubs()) {
        		
        		// Known value of stubs
            	int nstubsOnBranch = stubs == null ? 0 : stubs.getNStubsOnBranch(nodeNr);
            	
            	
            	
            	// Number of spikes is nstubs+1, unless the sibling is a sampled ancestor, in which case it is nstubs
            	Node node = treeInput.get().getNode(nodeNr);
				int spikeSum = getNSpikes(node, nstubsOnBranch);
				if (spikeSum == 0) {
					
					
					// Delta function
					double logprob = 0 ;
					if (spikeOfBranch != 0) {
						logprob = Double.NEGATIVE_INFINITY;
					}else {
						logprob = 0;
					}
					
					logP += logprob;
					
				}
				
				else {
            	
	            	double alphaBranch = shape * spikeSum; 
	            	gamma = new GammaDistributionImpl(alphaBranch, scale);
	            	logP += gamma.logDensity(spikeOfBranch);
            	
				}
        		
        	}else {
        		
        		// Unknown value - integrate across all possible values
        		Node node = treeInput.get().getNode(nodeNr);
        		double h0 = node.getHeight();
        		double h1 = node.isRoot() ? h0 : node.getParent().getHeight();
        		double mu = stubs.getMeanNumberOfStubs(h0, h1);
        		
        		//Log.warning("no est -> " + mu);
        		
        		
        		if (mu > 0) {
        			
        			double branchP = 0;
        			int k = 0;
        			double cumsum = 0;
        			while (cumsum < MAX_CUM_SUM) {
        				
        				// P(k observations) under a Poisson(mu)
        				double p = -mu + k*Math.log(mu);
        				for (int i = 2; i <= k; i ++) p += -Math.log(i);
        				double pReal = Math.exp(p);
        				
        				
        				cumsum += pReal;
        				
        				
        				// Number of spikes is nstubs+1, unless the sibling is a sampled ancestor, in which case it is nstubs
        				int spikeSum = getNSpikes(node, k);
        				if (spikeSum == 0) {
        					
        					
        					// Delta function
        					if (spikeOfBranch != 0) {
        						branchP += 0;
        					}else {
        						branchP += Math.exp(p);
        					}
        					
        				}
        				
        				else {
        				
	        				double alphaBranch = shape * spikeSum; 
	        				gamma = new GammaDistributionImpl(alphaBranch, scale);
	        				double gammaLogP = gamma.logDensity(spikeOfBranch);
	        				if (spikeOfBranch == 0|| gammaLogP == Double.NEGATIVE_INFINITY || Double.isNaN(gammaLogP)) {
	        					branchP += 0;
	        				}else {
	        					branchP += Math.exp(p + gammaLogP);
	        				}
        				
        				}
        				
        				
        				k++;
        				
        			}
        			
        			//Log.warning(spikeOfBranch + " k = " + (k-1) + " -> " + Math.log(branchP));
        			
        			
        			logP += Math.log(branchP);
        			
        		}
        		
        		else {
        			
        			
        			int spikeSum = getNSpikes(node, 0);
    				if (spikeSum == 0) {
	        			
	        			//Log.warning("SA " + spikeOfBranch);
	        			
	        			// Delta function
						if (spikeOfBranch != 0) {
							logP += Double.NEGATIVE_INFINITY;
						}else {
							logP += 0;
						}
						
    				}else {
    					
    					//Log.warning("not SA " + spikeOfBranch);
        			
    					gamma = new GammaDistributionImpl(shape, scale);
    					logP += gamma.logDensity(spikeOfBranch);
    					
    				}
        			
        		}
        		
        		
        		
        	}
        
        	
        }
        
        // Numerical issue
        if (logP == Double.POSITIVE_INFINITY) {
        	logP = Double.NEGATIVE_INFINITY;
        }
        
        if (logP == Double.NEGATIVE_INFINITY) {
        	//Log.warning("Ninf");
        }
        
        
        return logP;
    }
	
	
	/**
	 * Number of spikes on a branch is number of stubs plus 1
	 * Unless the parent is a sampled ancestor, in which case it is just the number of stubs
	 * @param node
	 * @param nstubs
	 * @return
	 */
	public static int getNSpikes(Node node, int nstubs) {
		
		if (node.isRoot()) {
			return nstubs + 1;
		}
		
		if (node.isDirectAncestor()) {
			return 0;
		}
		
		if (node.getParent().isFake()) {
			return nstubs;
		}
		
		return nstubs + 1;
		
		
		
	}
	
	@Override
	public List<String> getConditions() {
		List<String> conds = new ArrayList<>();
		//conds.add(meanInput.get().getID());
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
        double mean = 1; //meanInput.get().getValue();
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
        	
        	
        	Node node = treeInput.get().getNode(nodeNr);;
        	int nspikes = getNSpikes(node, nstubsOnBranch);
        	
        	if (nspikes == 0) {
        		double spikeOfBranch = 0;
        		spikesInput.get().setValue(nodeNr, spikeOfBranch);
        	}
        	
        	else {
        		
        		double alphaBranch = shape * nspikes; // One spike for the branch, and one per stub
            	gamma = new GammaDistributionImpl(alphaBranch, scale);
            	
            	try {
    				double spikeOfBranch = gamma.inverseCumulativeProbability(random.nextFloat());
    				spikesInput.get().setValue(nodeNr, spikeOfBranch);
    				
    			} catch (MathException e) {
    				e.printStackTrace();
    				throw new IllegalArgumentException("Unexpected error when sampling from Gamma(" + shape + ", " + scale + ")");
    			}
        	}
        	
        	
        	
        }
		
		
		
		
	}
	@Override
    protected boolean requiresRecalculation() {
        return super.requiresRecalculation() || 
        		InputUtil.isDirty(stubsInput) || 
        		InputUtil.isDirty(spikesInput) || 
        		InputUtil.isDirty(shapeInput) || 
        		InputUtil.isDirty(meanInput);
    }

	
	
	/**
	 * Calculate cumulative probabilities of sampling a stub, conditional on the gamma distribution and tree prior (theta)
	 * p(nstubs | spike size, theta) = p(spike size | nstubs, theta) x p (nstubs, theta) / p (spike size | theta)
	 * @param h0
	 * @param h1
	 * @return
	 */
	public double[] getCumulativeProbs(double mu, int nodeNr) {
		
		
		List<Double> probs = new ArrayList<>();
		
		
		
		// Shape and scale of gamma 
		double shape = shapeInput.get().getValue();
		double mean = 1; //meanInput.get().getValue();
		if (shape <= 0 || mean <= 0) {
			throw new IllegalArgumentException("Cannot sample spikes because shape or mean are non-positive " + shape + "  " + mean);
		}
		double scale = mean / shape;
		
		// Spike size
		double spikeOfBranch = spikesInput.get().getValue(nodeNr);
		//if (meanGamma1Input.get()) {
    		//scale = 1 / shape;
    	//}

		
		
		int k = 0;
		double poissonCumSum = 0;
		
		double branchPSum = 0;
		while (poissonCumSum < MAX_CUM_SUM) {
			
			double branchP = 0;
			
			// P(k observations) under a Poisson(mu)
			double p = -mu + k*Math.log(mu);
			for (int i = 2; i <= k; i ++) p += -Math.log(i);
			double pReal = Math.exp(p);
			
			
			
			// Integrate across all possible values in poisson distribution
			double alphaBranch = shape * (k + 1); // One spike for the branch, and one per stub
			
			
			if (indicatorInput.get() != null && indicatorInput.get().getValue()) {
			
				gamma = new GammaDistributionImpl(alphaBranch, scale);
				double gammaLogP = gamma.logDensity(spikeOfBranch);
				if (gammaLogP == Double.NEGATIVE_INFINITY || Double.isNaN(gammaLogP)) {
					if (poissonCumSum > 0) break;
					branchP += 0;
				}else {
					branchP += Math.exp(p + gammaLogP);
				}
				
				//Log.warning(mu + ": k=" + k + " -> " + gammaLogP);
				//Log.warning("" + Math.exp(p + gammaLogP));
				
			}else {
				branchP += pReal;
			}
			
			
			poissonCumSum += pReal;
			branchPSum += branchP;
			probs.add(branchP);
			
			k++;
		}
		
		
		// Normalise to sum to 1
		double[] array = new double[probs.size()];
		for(int i = 0; i < probs.size(); i++) {
			array[i] = probs.get(i) / branchPSum; 
		}
		
		
		// Convert into cumulative sum
		double cumsum = 0;
		for(int i = 0; i < probs.size(); i++) {
			double p = array[i];
			array[i] = p + cumsum;
			cumsum += p;
			//Log.warning("P(K<=" + i + ") = " + array[i]);
		}
		return array;
	}
	

}
