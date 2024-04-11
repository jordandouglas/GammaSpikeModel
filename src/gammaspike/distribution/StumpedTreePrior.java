package gammaspike.distribution;


import java.io.PrintStream;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.core.Input.Validate;
import beast.base.evolution.speciation.SpeciesTreeDistribution;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.evolution.tree.TreeIntervals;
import beast.base.evolution.tree.TreeUtils;
import beast.base.inference.distribution.Poisson;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import gammaspike.tree.Stubs;


// TODO: does not give correct posterior when origin is used
@Description("Prior distribution on a stumped tree")
public class StumpedTreePrior extends SpeciesTreeDistribution {

	final public Input<RealParameter> lambdaInput = new Input<>("lambda", "birth rate lambda", Validate.REQUIRED);
	//final public Input<RealParameter> muInput = new Input<>("mu", "death rate mu", Validate.REQUIRED);
	final public Input<RealParameter> r0Input = new Input<>("r0", "ratio between birth and death rate", Input.Validate.REQUIRED);
	final public Input<RealParameter> originInput = new Input<>("origin", "length of origin branch", Validate.OPTIONAL);
	final public Input<Stubs> stubsInput = new Input<>("stubs", "the stubs of this tree", Input.Validate.REQUIRED);
	final public Input<Integer> maxNrStubsInput = new Input<>("maxNrStubs", "max number of stubs", -1);
	final public Input<Boolean> ignoreTreePriorInput = new Input<>("ignoreTreePrior", "ignore tree prior (debugging)", false);
	
	final public Input<RealParameter> pInput = new Input<>("p", "the probability of a single-mutant. if included, then stubs correspond to bursts, and"
			+ " the rate of stub appearance is multiplied by 1-p/2", Input.Validate.OPTIONAL);
	final public Input<Boolean> perBranchSpikeInput = new Input<>("perBranchSpike", "if true, then the indicator will disregard labels", false);

    
	final double EPSILON = 1e-8; // Numerical error 
	
	
	boolean initialising = true;
	
	@Override
    public void initAndValidate() {
        super.initAndValidate();
        
        TreeInterface tree = treeInput.get();
        if (tree == null) {
            tree = treeIntervalsInput.get().treeInput.get();
        }
        if (!TreeUtils.isUltrametric(tree)) {
            Log.warning.println("WARNING: This model (tree prior) cannot handle dated tips. " +
                    "Please select a tree prior that can, otherwise results may be invalid.");
        }
        
        
        
        
    }
	
	
	@Override
    public double calculateTreeLogLikelihood(final TreeInterface treeInterface) {
        
		
		
		// Get variables
		double p = 0;
		Tree tree = (Tree) treeInterface;
		double lambda = lambdaInput.get().getValue();
		double mu = lambda / r0Input.get().getValue();
		Double originLength = originInput.get() == null ? null : originInput.get().getValue();
		Stubs stubs = stubsInput.get();
		
		//TreeIntervals intervals = new TreeIntervals(tree);
		
		
		if (maxNrStubsInput.get() > -1 && stubs.getStubDimension() > maxNrStubsInput.get()) {
			if (initialising) Log.warning("Cannot initialise because there are too many stubs");
			initialising = false;
			return Double.NEGATIVE_INFINITY;
		}
		
		if (ignoreTreePriorInput.get()) {
			return 0;
		}
		
		
		// Condition checker
		if (originLength != null && originLength < 0) {
			if (initialising) Log.warning("Cannot initialise because origin length is less than 0");
			initialising = false;
			return Double.NEGATIVE_INFINITY;
		}
		if (lambda <= mu) {
			if (initialising) Log.warning("Cannot initialise because lambda is less than mu");
			initialising = false;
			return Double.NEGATIVE_INFINITY;
		}
		
		double blueLogP = 0;
		double redLogP = 0;
		
		
		
		
		
		// Blue tree speciation density, conditional on survival of both children
		blueLogP = -this.getBlueIntegral();
		for (int branchNr = 0; branchNr < tree.getNodeCount(); branchNr++) {
			if (tree.getNode(branchNr).isLeaf()) continue;
			double height = tree.getNode(branchNr).getHeight();
			double logRate = this.getBlueTreeLogBirthRate(height, lambda, mu);
			blueLogP += logRate;
		}
		
		
		//if (true) return blueLogP;
		
		
		
		
		// Red tree stub density
		for (int branchNr = 0; branchNr < tree.getNodeCount(); branchNr++) {
			redLogP += getLogPForBranch(branchNr);
		}
		
	
		
		if (initialising && blueLogP == Double.NEGATIVE_INFINITY) {
			Log.warning("Cannot initialise stumped tree prior because blue p is negative infinity");
		}
		
		else if (initialising && redLogP == Double.NEGATIVE_INFINITY) {
			Log.warning("Cannot initialise stumped tree prior because red p is negative infinity");
		}
		
		
		p = blueLogP + redLogP;
		initialising = false;
		return p;
		
		
		
    }
	
	
	/**
	 * Red tree 
	 * @param branchNr
	 * @return
	 */
	public double getLogPForBranch(int branchNr) {
		

		// Get variables
		double p = 0;
		Tree tree = (Tree) treeInput.get();
		double lambda = lambdaInput.get().getValue();
		double mu = lambda / r0Input.get().getValue();
		Double originLength = originInput.get() == null ? null : originInput.get().getValue();
		Stubs stubs = stubsInput.get();
		double T = tree.getRoot().getHeight() + originLength;
		
		
		Node node = tree.getNode(branchNr);
		
		// Reverse times
		double h0 = node.isRoot() ? T : node.getParent().getHeight();
		double h1 = node.getHeight();
		
		// Stubs on this branch
		int m = 0;
		for (int stubNr = 0; stubNr < stubs.getStubDimension(); stubNr++) {
			if (!stubs.includeStub(stubNr)) continue;
			if (stubs.getBranches().getValue(stubNr) == branchNr) {
				double stubHeight = stubs.getAbsoluteTimeOfStub(stubNr);
				//double stubHeight = stubs.getRelativeTimeOfStub(stubNr);
				
				if (stubHeight < node.getHeight() || stubHeight > node.getParent().getHeight()) {
					return Double.NEGATIVE_INFINITY;
				}
				double logRate = this.getRedTreeBirthLogRate(stubHeight, lambda, mu);
				p += logRate;
				m++;
			}
		}
		
		
		// Integral and combinatorial term across the branch
		double g = getGForInterval(h1, h0);
		p += -g;
		for (int i = 2; i <= m; i ++) p += -Math.log(i);
		
		return p;
		
	}
	

	
	
	public double getGForInterval(double h1, double h0) {
		
		Double originLength = originInput.get() == null ? null : originInput.get().getValue();
		Tree tree = (Tree) treeInput.get();
		double lambda = lambdaInput.get().getValue();
		double mu = lambda / r0Input.get().getValue();
		double T = tree.getRoot().getHeight() + originLength;
		
		return getRedRateIntegral(h1, lambda, mu, T) - getRedRateIntegral(h0, lambda, mu, T);
	}

//
//	/**
//	 * How many lineages cross this time?
//	 * @param tree
//	 * @param height
//	 * @return
//	 */
//	private int getNrBranchesAtHeight(Tree tree, double height, double originLength) {
//		
//		int b = 0;
//		for (Node node : tree.getNodesAsArray()) {
//			
//			double h1 = node.getHeight();
//			double h2 = node.isRoot() ? (node.getHeight() + originLength) : node.getParent().getHeight();
//			if (height >= h1 & height < h2) b++;
//			
//		}
//		
//		
//		//Log.warning("at height " + height + " there are " + b);
//		
//		return b;
//		
//	}
//	
//	
	
	
	/**
	 * Multiply stub rate by this term
	 * @return
	 */
	public double getStubMultiplier() {
		if (pInput.get() == null) return 1.0;
		if (pInput.get().getDimension() == 1 && !perBranchSpikeInput.get()) {
			return 1-pInput.get().getValue()/2;
		}else if (pInput.get().getDimension() == 1 && perBranchSpikeInput.get()) {
			return pInput.get().getValue();
		}else {
			double p1 = pInput.get().getValue(0); // Single mutant
			double p2 = pInput.get().getValue(1); // Double mutant
			double p3 = pInput.get().getValue(2); // No mutant
			return p1/2 + p2;
		}
		
	}

	
	
	/**
	 * Red rate integrated from top to bottom of tree
	 * @return
	 */
	public double getG() {
		
		Tree tree = (Tree) treeInput.get();
		double lambda = lambdaInput.get().getValue();
		double mu = lambda / r0Input.get().getValue();
		Double originLength = originInput.get() == null ? null : originInput.get().getValue();
		double T = tree.getRoot().getHeight() + originLength;
		
		
		double Q = Math.log(lambda-mu) -Math.log(lambda*Math.exp((lambda-mu)*T) - mu);
		double g = 2*Q + 2*lambda*T;
		return g;
		
	}
	
	
	
	/**
	 * Red rate integrated along a branch
	 * @return
	 */
	public double getG(int branchNr) {
		
		
		Tree tree = (Tree) treeInput.get();
		double lambda = lambdaInput.get().getValue();
		double mu = lambda / r0Input.get().getValue();
		Double originLength = originInput.get() == null ? null : originInput.get().getValue();
		double T = tree.getRoot().getHeight() + originLength;
		

			
		Node node = tree.getNode(branchNr);
		
		// Reverse times
		double h0 = node.isRoot() ? (node.getHeight() + originLength) : node.getParent().getHeight();
		double h1 = node.getHeight();
		
		// Integral term across the branch
		double g = getRedRateIntegral(h1, lambda, mu, T) - getRedRateIntegral(h0, lambda, mu, T);
		return g;
		
	}
	
	
	private double getGSum() {
		
		double gSum = 0;
		
		Tree tree = (Tree) treeInput.get();
		double lambda = lambdaInput.get().getValue();
		double mu = lambda / r0Input.get().getValue();
		Double originLength = originInput.get() == null ? null : originInput.get().getValue();
		double T = tree.getRoot().getHeight() + originLength;
		
		// For each branch
		for (int branchNr = 0; branchNr < tree.getNodeCount(); branchNr++) {
			
			
			Node node = tree.getNode(branchNr);
			
			// Reverse times
			double h0 = node.isRoot() ? (node.getHeight() + originLength) : node.getParent().getHeight();
			double h1 = node.getHeight();
			
			// Integral term across the branch
			double g = getRedRateIntegral(h1, lambda, mu, T) - getRedRateIntegral(h0, lambda, mu, T);
			gSum += g;
			
		}
		
		return gSum;
	}
	
	/**
	 * Red term integrated across the whole tree - aka expected number of stubs
	 * @param tree
	 * @param lambda
	 * @param mu
	 * @return
	 */
	private double getH() {
		
		Tree tree = (Tree) treeInput.get();
		double lambda = lambdaInput.get().getValue();
		double mu = lambda / r0Input.get().getValue();
		Double originLength = originInput.get() == null ? 0 : originInput.get().getValue();
		double T = tree.getRoot().getHeight();
		
		
		double g = 0;
		for (Node node : tree.getNodesAsArray()) {
			
			
			// Reverse times
			double h0 = node.isRoot() ? (node.getHeight() + originLength) : node.getParent().getHeight();
			double h1 = node.getHeight();
			double integralBranch = getRedRateIntegral(h1, lambda, mu, T) - getRedRateIntegral(h0, lambda, mu, T);
			g += integralBranch;
			
		}
		
		
		return g;
		
	}
	
	
	/**
	 * Blue rate integrated across the whole tree
	 * @param tree
	 * @param lambda
	 * @param mu
	 * @return
	 */
	private double getBlueIntegral() {
		
		Tree tree = (Tree) treeInput.get();
		double lambda = lambdaInput.get().getValue();
		double mu = lambda / r0Input.get().getValue();
		Double originLength = originInput.get() == null ? 0 : originInput.get().getValue();
		double T = tree.getRoot().getHeight();
		
		
		double blueIntegral = 0;
		for (Node node : tree.getNodesAsArray()) {
			
			
			// Reverse times
			double h0 = node.isRoot() ? (node.getHeight() + originLength) : node.getParent().getHeight();
			double h1 = node.getHeight();
			double integralBranch = getBlueRateIntegral(h1, lambda, mu, T) - getBlueRateIntegral(h0, lambda, mu, T);
			blueIntegral += integralBranch;
			
		}
		
		
		return blueIntegral;
		
	}
	
	
	private double getBlueRateIntegral(double height, double lambda, double mu, double rootHeight) {
		double finalTime = rootHeight;
		double t = rootHeight - height;
		//Log.warning(t +  " " + finalTime + " " + lambda + " " + mu);
		double tlm = t*(lambda - mu);
		double b = Math.log(lambda * Math.exp(lambda * finalTime) - mu*Math.exp(tlm + mu*finalTime));
		return tlm - b;
	}
	
	
	private double getRedRateIntegral(double height, double lambda, double mu, double finalTime) {
		double s = finalTime - height;
		double b = Math.log(lambda*Math.exp( (lambda-mu)*(finalTime - s)) - mu);
		return getStubMultiplier()*2*(b + lambda*s);
	}
	
	
	private double getBlueTreeLogBirthRate(double timeRemaining, double lambda, double mu) {
		double logqt = this.getLogQt(timeRemaining, lambda, mu);
		return Math.log(lambda) + Math.log(1 - Math.exp(logqt));
	}
	
	
	
	private double getRedTreeBirthLogRate(double timeRemaining, double lambda, double mu) {
		double logqt = this.getLogQt(timeRemaining, lambda, mu);
		return Math.log(getStubMultiplier()) + Math.log(2) + Math.log(lambda) +  logqt;
	}
	
	
	
	/**
	 * Get qt
	 * @param time
	 * @param lambda
	 * @param mu
	 * @return
	 */
	private double getLogQt(double time, double lambda, double mu) {
		double exp = Math.exp((lambda - mu)*time);
		double top = (exp - 1);
		double bottom = lambda*exp - mu;
		double logqt = Math.log(mu) + Math.log(top) - Math.log(bottom);
		return logqt;
	}
	
	
	/*
	private double getRedTreeBirthRate(double timeRemaining, double lambda, double mu) {
		double qt = this.getQt(timeRemaining, lambda, mu);
		return 2 * lambda * qt;
		//return 2 * lambda * qt * (1 - qt);
	}
	

	private double getQt(double time, double lambda, double mu) {
		
		
		double exp = Math.exp((lambda - mu)*time);
		double top = mu * (exp - 1);
		double bottom = lambda*exp - mu;
		double qt = top / bottom;
		
		
		// Round to 0 or 1 if error is small
		if (qt < 0 && qt > -EPSILON) qt = 0;
		else if (qt > 1 && qt < 1+EPSILON) qt = 1;
		
		if (qt < 0 | qt > 1) {
			Log.warning("qt is not a probability " + qt);
		}
		//Log.warning("qt = " + qt);
		return qt;
		
	}
	*/
	
	@Override
    protected boolean requiresRecalculation() {
		if (originInput.get() != null && originInput.get().somethingIsDirty()) {
			return true;
		}
		if (pInput.get() != null && pInput.get().somethingIsDirty()) {
			return true;
		}
        return super.requiresRecalculation() || lambdaInput.get().somethingIsDirty() || r0Input.get().somethingIsDirty() || InputUtil.isDirty(stubsInput);
    }

	
    /**
     * Loggable interface implementation follows *
     */
    @Override
    public void init(final PrintStream out) {
       // out.print(getID() + "\t" + "g0\tg1\tg2\tg3\t"+ "stubGsum\t"+ "stubH\t");
        out.print(getID() + "\t" + "stubGsum\t"+ "stubH\t");
    }

    @Override
    public void log(final long sample, final PrintStream out) {
    	//out.print(getCurrentLogP() + "\t" + this.getG(0) + "\t" + this.getG(1) + "\t" + this.getG(2) + "\t" + this.getG(3) + "\t" + this.getGSum()  +"\t"+ this.getH() + "\t");
        out.print(getCurrentLogP() + "\t" + this.getGSum()  +"\t"+ this.getH() + "\t");
    }


	public double getTotalTime() {
		Tree tree = (Tree) treeInput.get();
		return tree.getRoot().getHeight();
	}


	public Tree getTree() {
		return (Tree) treeInput.get();
	}

	
	

	
	// CDF of this time (not height) in this time interval
	// y is the time (not height)
	// t0 is the time of the start of the interval
	// t1 is the time of the end
	// T is the total time (aka root height)
	public double getCDF(double y, double t0, double t1) {
		
		
		double T = this.getTotalTime();
		double lambda = lambdaInput.get().getValue();
		double mu = lambda / r0Input.get().getValue();
		
		
		if (y <= t0) return 0;
		if (y >= t1) return 1;
		
		double f0 = getRedRateIntegral(T-t0, lambda, mu, T);
		double f1 = getRedRateIntegral(T-t1, lambda, mu, T);
		double fy = getRedRateIntegral(T-y, lambda, mu, T);
		
		double p = (fy - f0) / (f1 - f0);
		
		
		if (p < -EPSILON || p > 1 + EPSILON) {
			Log.warning("cdf p = " + p + " for y = " + y + " " + f0 + " " + f1 + " " + fy);
		}
		
		if (p < 0) p = 0;
		if (p > 1) p = 1;
		
		return p;
		
		
	}
	
	

	// Sample a stub time along [0,T] interval using a numeric approximation of the icdf (bisection method)
	// t0 is the start
	// t1 is the end
	public double sampleTime(double t0, double t1) {
		
		
		
		double T = this.getTotalTime();
		
		
		if (t0 > t1 || t0 > T || t1 > T) throw new IllegalArgumentException("Dev error 21414: invalid t0 or t1");
		
		
		double a = t0;
		double b = t1;
		double y = (a+b)/2;
		double u = Randomizer.nextDouble();
		double cdf = this.getCDF(y, t0, t1);
		
		// Keep changing y until cdf(y) ~= u
		int ninter = 0;
		while (Math.abs(cdf - u) > EPSILON) {
			
			
			if (cdf < u) {
				
				// Increase y
				a = y;
				
			}else {
				
				// Decrease y
				b = y;
				
			}
			
			y = (a+b)/2;
			
			cdf = this.getCDF(y, t0, t1);
			
			
			ninter++;
			if (ninter > 10000) {
				Log.warning("Warning: cannot solve for u=" + u + ", cdf(y=" + y + ")=" + cdf);
				break;
			}
		}
		
		
		//Log.warning("Warning: successfully solved for u=" + u + ", cdf(y=" + y + ")=" + cdf);
		
		
		// Warning: this is a time not a height!
		return y;
	}


	
	
}


