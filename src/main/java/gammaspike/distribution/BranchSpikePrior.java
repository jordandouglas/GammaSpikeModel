package gammaspike.distribution;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import org.apache.commons.statistics.distribution.GammaDistribution;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.spec.domain.NonNegativeInt;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.inference.distribution.TensorDistribution;
import beast.base.spec.inference.parameter.BoolScalarParam;
import beast.base.spec.inference.parameter.IntVectorParam;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.type.RealVector;
import gammaspike.tree.Stubs;
import beast.base.util.Randomizer;





@Description("A sum of gamma distributions, one for each spike on a branch")
public class BranchSpikePrior extends TensorDistribution<RealVector<NonNegativeReal>, Double> {
	
	final public Input<Stubs> stubsInput = new Input<>("stubs", "stubs of the tree", Input.Validate.OPTIONAL);
	final public Input<IntVectorParam<? extends NonNegativeInt>> nstubsInput = new Input<>("nstubs", "num stubs per branch", Input.Validate.OPTIONAL);
	final public Input<RealScalarParam<? extends PositiveReal>> shapeInput = new Input<>("shape", "shape parameter for the gamma distribution of each spike.", Input.Validate.REQUIRED);
	final public Input<RealScalarParam<? extends PositiveReal>> meanInput = new Input<>("mean", "mean (=shape*scale) parameter for the gamma distribution of each spike.", Input.Validate.OPTIONAL); 
	
	final public Input<Tree> treeInput = new Input<>("tree", "tree required for setting the spike dimension (if direct sampling)", Input.Validate.OPTIONAL); 
	final public Input<BoolScalarParam> indicatorInput = new Input<>("indicator", "burst size is 0 if this is false", Input.Validate.OPTIONAL);

	
	// If there are too many stubs on a branch (e.g., during mixing) then the gamma distribution shape is large, which causes instabilities
	final double MAX_CUM_SUM = 0.999;
	
	
	
	@Override
    public void initAndValidate() {
		 
		if (stubsInput.get() != null) {
			stubsInput.get().setBranchSpikePrior(this);
		}
		 
	}
	
	@Override
	public double calculateLogP() {
		
		
		RealVector<NonNegativeReal> branchRates = paramInput.get();
		logP = calcLogP(branchRates.getElements());
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
		//args.add(paramInput.get().getID());
		return args;
	}
	
	
	@Override
	public List<Double> sample() {
		

		if (treeInput.get() == null) {
			throw new IllegalArgumentException("Please specify the tree");
		}
		
		
		List<Double> vals = new ArrayList<>();

		Tree tree = (Tree) treeInput.get();
		int dimension = tree.getNodeCount();
		
		
	    // Check shape and scale are positive
        double shape = shapeInput.get().get();
        double mean = 1; //meanInput.get().getValue();
        if (shape <= 0 || mean <= 0) {
        	throw new IllegalArgumentException("Cannot sample spikes because shape or mean are non-positive " + shape + "  " + mean);
    	}
        double scale = mean / shape;
        
        // Calculate density of the total spike size of each branch, assuming that each node or stub has an iid spike drawn from a Gamma(alpha, beta)
        // This approach integrates across all stub spike sizes, so we don't need to estimate them individually
        Stubs stubs = stubsInput.get();
        for (int nodeNr = 0; nodeNr < dimension; nodeNr ++) {
        	
        	int nstubsOnBranch = 0;
        	if (nodeNr < nstubsInput.get().size()) {
        		nstubsOnBranch = stubs == null ? nstubsInput.get().get(nodeNr) : stubs.getNStubsOnBranch(nodeNr);
        	}
        	
        	Node node = treeInput.get().getNode(nodeNr);;
        	int nspikes = getNSpikes(node, nstubsOnBranch);
        	
        	if (nspikes == 0) {
        		double spikeOfBranch = 0;
        		vals.add(spikeOfBranch);
        	} else {
        		double alphaBranch = shape * nspikes; // One spike for the branch, and one per stub
        		
        		GammaDistribution gamma = GammaDistribution.of(alphaBranch, scale);
            	try {
    				double spikeOfBranch = gamma.inverseCumulativeProbability(Randomizer.nextFloat());
    				vals.add(spikeOfBranch);
    			} catch (Exception e) {
    				e.printStackTrace();
    				throw new IllegalArgumentException("Unexpected error when sampling from Gamma(" + shape + ", " + scale + ")");
    			}
        	}

        }
        
        return vals;
		
		
	}
	
	
//
//	@Override
//	// Sample a new "spike" value for every node (or branch) in a tree
//	public void sample(State state, Random random) {
//		
//		if (treeInput.get() == null) {
//			throw new IllegalArgumentException("Please specify the tree");
//		}
//		
//		if (sampledFlag) return;
//		sampledFlag = true;
//
//		// Cause conditional parameters to be sampled
//		sampleConditions(state, random);
//
//		Tree tree = (Tree) treeInput.get();
//		int dimension = tree.getNodeCount();
//		spikesInput.get().setDimension(dimension);
//		//spikesInput.get().setValue(null);
//		
//	    // Check shape and scale are positive
//        double shape = shapeInput.get().get();
//        double mean = 1; //meanInput.get().getValue();
//        if (shape <= 0 || mean <= 0) {
//        	throw new IllegalArgumentException("Cannot sample spikes because shape or mean are non-positive " + shape + "  " + mean);
//    	}
//        double scale = mean / shape;
//        
//        // Calculate density of the total spike size of each branch, assuming that each node or stub has an iid spike drawn from a Gamma(alpha, beta)
//        // This approach integrates across all stub spike sizes, so we don't need to estimate them individually
//        Stubs stubs = stubsInput.get();
//        for (int nodeNr = 0; nodeNr < spikesInput.get().size(); nodeNr ++) {
//        	
//        	int nstubsOnBranch = 0;
//        	if (nodeNr < nstubsInput.get().size()) {
//        		nstubsOnBranch = stubs == null ? nstubsInput.get().get(nodeNr) : stubs.getNStubsOnBranch(nodeNr);
//        	}
//        	
//        	Node node = treeInput.get().getNode(nodeNr);;
//        	int nspikes = getNSpikes(node, nstubsOnBranch);
//        	
//        	if (nspikes == 0) {
//        		double spikeOfBranch = 0;
//        		spikesInput.get().set(nodeNr, spikeOfBranch);
//        	} else {
//        		double alphaBranch = shape * nspikes; // One spike for the branch, and one per stub
//        		
//        		GammaDistribution gamma = GammaDistribution.of(alphaBranch, scale);
//            	try {
//    				double spikeOfBranch = gamma.inverseCumulativeProbability(random.nextFloat());
//    				spikesInput.get().set(nodeNr, spikeOfBranch);
//    			} catch (Exception e) {
//    				e.printStackTrace();
//    				throw new IllegalArgumentException("Unexpected error when sampling from Gamma(" + shape + ", " + scale + ")");
//    			}
//        	}
//
//        }
//
//	}


	@Override
    protected boolean requiresRecalculation() {
		return true;
//        return super.requiresRecalculation() ||
//        		InputUtil.isDirty(stubsInput) ||
//        		InputUtil.isDirty(spikesInput) ||
//        		InputUtil.isDirty(shapeInput) ||
//        		InputUtil.isDirty(meanInput);
    }

	
	
	/**
	 * Calculate cumulative probabilities of sampling a stub, conditional on the gamma distribution and tree prior (theta)
	 * p(nstubs | spike size, theta) = p(spike size | nstubs, theta) x p (nstubs | theta) / p (spike size | theta)
	 * @param mu (Poisson distribution mean)
	 * @param nodeNr
	 * @return
	 */
	// Calculates the cumulative posterior probability distribution for the number of "stubs" on a specific branch (nodeNr)
	// "Given the spikeOfBranch value we observed, what is the probability that this branch has 0 stubs, ≤1 stub, ≤2 stubs, etc.?"
	// The method is used by sampleNStubsOnBranch (in Stubs.java) to sample the number of stubs for a branch
	public double[] getCumulativeProbs(double mu, int nodeNr) {
		
		List<Double> probs = new ArrayList<>();
		
		// Shape and scale of gamma
		double shape = shapeInput.get().get();
		double mean = 1;
		if (shape <= 0) {
			throw new IllegalArgumentException("Cannot sample spikes because shape or mean are non-positive " + shape + "  " + mean);
		}
		double scale = mean / shape;
		
		// Non-weighted spike size (spike mean not taken into account)
		double spikeOfBranch = paramInput.get().get(nodeNr);

		int k = 0;
		double poissonCumSum = 0;
		
		double branchPSum = 0;
		while (poissonCumSum < MAX_CUM_SUM) {
			
			double branchP = 0;
			
			// P(k observations) under a Poisson(mu)
			double p = -mu + k * Math.log(mu);
			for (int i = 2; i <= k; i ++) p += -Math.log(i);
			double pReal = Math.exp(p);
			
			// Integrate across all possible values in poisson distribution
			double alphaBranch = shape * (k + 1); // One spike for the branch, and one per stub
			
			// If the use-spike indicator is true
			if (indicatorInput.get() != null && indicatorInput.get().get()) {
			
				GammaDistribution gamma = GammaDistribution.of(alphaBranch, scale);
				//gamma = new Gamma();
            	//gamma.initByName("alpha", alphaBranch, "theta", scale);
				double gammaLogP = gamma.logDensity(spikeOfBranch);
				if (gammaLogP == Double.NEGATIVE_INFINITY || Double.isNaN(gammaLogP)) {
					if (poissonCumSum > 0) break;
					branchP += 0;
				} else {
					branchP += Math.exp(p + gammaLogP);
				}
				
			// If the use-spike indicator is false or not provided
			// It ignores the spikeOfBranch value completely and just uses the prior: branchP = P(K=k).
			} else {
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

	@Override
	public void refresh() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public Double getLowerBoundOfParameter() {
		return 0.0;
	}

	@Override
	public Double getUpperBoundOfParameter() {
		return Double.POSITIVE_INFINITY;
	}

	@Override
	protected double calcLogP(Double... value) {
		return this.calcLogP(Arrays.asList(value));
	}
	
	
    private double calcLogP(List<Double> branchRates) {
    	
    	
    	double logp = 0;
    	
    	
        // Check spike shape and scale are positive
        double shape = shapeInput.get().get();
        double mean = 1; // meanInput.get().getValue();
        if (shape <= 0 || mean <= 0) {
    		return Double.NEGATIVE_INFINITY;
    	}
        double scale = mean / shape;
        
        // Calculate density of the total spike size of each branch, assuming that each node or stub has an iid spike drawn from a Gamma(alpha, beta)
        // This approach integrates across all stub spike sizes, so we don't need to estimate them individually
        Stubs stubs = stubsInput.get();
        for (int nodeNr = 0; nodeNr < paramInput.get().size(); nodeNr ++) {
        	
        	double spikeOfBranch = paramInput.get().get(nodeNr);
        	if (spikeOfBranch < 0) {
        		return Double.NEGATIVE_INFINITY;
        	}

        	if (stubs == null || stubs.estimateStubs()) { // Either stubsInput is null or stubs are being estimated
        		
        		// Known value of stubs
				// If stubsInput is null, nstubsOnBranch is zero; otherwise get the number of stubs on branch nodeNr
            	int nstubsOnBranch = stubs == null ? 0 : stubs.getNStubsOnBranch(nodeNr);
            	
            	// Number of spikes is nstubs + 1, unless the sibling is a sampled ancestor, in which case it is nstubs
            	Node node = treeInput.get().getNode(nodeNr);
				int spikeSum = getNSpikes(node, nstubsOnBranch);
				if (spikeSum == 0) {
					// Delta function
					double logprob = 0 ;
					if (spikeOfBranch != 0) {
						logprob = Double.NEGATIVE_INFINITY;
					} else {
						logprob = 0;
					}
					logp += logprob;
				}
				else {
	            	double alphaBranch = shape * spikeSum;
	            	GammaDistribution gamma = GammaDistribution.of(alphaBranch, scale);
	            	logp += gamma.logDensity(spikeOfBranch);
				}
        		
        	} else { // Integrating over stubs (Stub-free inference)
        		
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
        				for (int i = 2; i <= k; i ++) p += -Math.log(i); // Integrating over all possible values
        				double pReal = Math.exp(p);

        				cumsum += pReal;
        				
        				// Number of spikes is nstubs + 1, unless the sibling is a sampled ancestor, in which case it is nstubs
        				int spikeSum = getNSpikes(node, k);
        				if (spikeSum == 0) {
        					// Delta function
        					if (spikeOfBranch != 0) {
        						branchP += 0;
        					} else {
        						branchP += Math.exp(p);
        					}
        				}
        				else {
	        				double alphaBranch = shape * spikeSum;
	        				GammaDistribution gamma = GammaDistribution.of(alphaBranch, scale);
	        				double gammaLogP = gamma.logDensity(spikeOfBranch);
	        				if (spikeOfBranch == 0|| gammaLogP == Double.NEGATIVE_INFINITY || Double.isNaN(gammaLogP)) {
	        					branchP += 0;
	        				} else {
	        					branchP += Math.exp(p + gammaLogP);
	        				}
        				}

        				k++;
        				
        			}
        			
        			logp += Math.log(branchP);
        			
        		}
        		
        		else {
        			
        			int spikeSum = getNSpikes(node, 0);
    				if (spikeSum == 0) {
	        			// Delta function
						if (spikeOfBranch != 0) {
							logp += Double.NEGATIVE_INFINITY;
						} else {
							logp += 0;
						}
    				} else {
    					GammaDistribution gamma = GammaDistribution.of(shape, scale);
    					logp += gamma.logDensity(spikeOfBranch);
    				}
        			
        		}

        	}

        }
        
        // Numerical issue
        if (logp == Double.POSITIVE_INFINITY) {
        	logp = Double.NEGATIVE_INFINITY;
        }

        
        return logp;
    	
    	
    	
    }
	

}
