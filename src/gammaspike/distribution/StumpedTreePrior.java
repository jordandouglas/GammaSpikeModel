package gammaspike.distribution;


import java.io.PrintStream;

import beast.base.core.Citation;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.core.Input.Validate;
import beast.base.evolution.speciation.SpeciesTreeDistribution;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.evolution.tree.TreeUtils;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import gammaspike.tree.Stubs;



@Description("Prior distribution on a stumped tree")
@Citation(value =
"Douglas, J., Bouckaert, R., Harris, S.C., Carter Jr, C.W., Wills, P.R. (2024) Evolution is coupled with branching across many granularities of life. bioRxiv 2024.09.08.611933", DOI = "https://doi.org/10.1101/2024.09.08.611933",
year = 2024, firstAuthorSurname = "Douglas")
public class StumpedTreePrior extends SpeciesTreeDistribution implements StubExpectation {

	final public Input<RealParameter> lambdaInput = new Input<>("lambda", "birth rate lambda", Validate.REQUIRED);
	//final public Input<RealParameter> muInput = new Input<>("mu", "death rate mu", Validate.REQUIRED);
	final public Input<RealParameter> r0Input = new Input<>("r0", "ratio between birth and death rate");
	final public Input<RealParameter> turnoverInput = new Input<>("turnover", "inverse of r0", Input.Validate.XOR, r0Input);
	
	final public Input<RealParameter> samplingProportionInput = new Input<>("samplingProportion", "sampling rate psi divided by (psi + mu)", Input.Validate.OPTIONAL);
	//final public Input<RealParameter> rhoInput = new Input<>("rho", "extant sampling probability rho (default 1)", Input.Validate.OPTIONAL);
	
	
	final public Input<Stubs> stubsInput = new Input<>("stubs", "the stubs of this tree", Input.Validate.OPTIONAL);
	final public Input<Boolean> ignoreTreePriorInput = new Input<>("ignoreTreePrior", "ignore tree prior (debugging)", false);
	final public Input<Boolean> ignoreStubPriorInput = new Input<>("ignoreStubPrior", "ignore tree prior (debugging)", false);
	
	final public Input<RealParameter> pInput = new Input<>("p", "the probability of a single-mutant. if included, then stubs correspond to bursts, and"
			+ " the rate of stub appearance is multiplied by 1-p/2", Input.Validate.OPTIONAL);
	final public Input<Boolean> perBranchSpikeInput = new Input<>("perBranchSpike", "if true, then the indicator will disregard labels", false);

    
	final double HEIGHT_THRESHOLD_OF_LEAF = 1e-10;
	int nExtantTaxa; // n
	
	boolean initialising = true;
	
	
	@Override
	public boolean canHandleTipDates() {
		return samplingProportionInput.get() != null;
	}
	
	@Override
    public void initAndValidate() {
        super.initAndValidate();
        
        TreeInterface tree = treeInput.get();
        if (tree == null) {
            tree = treeIntervalsInput.get().treeInput.get();
        }
       
        

        // THis dist is conditional on ntaxa, so ensure the number does not change
        nExtantTaxa = -1;
        
        
        
        // Ensure valid initial state
        final int nattempts = 1000;
		double lambda = lambdaInput.get().getValue();
		double mu = this.getMu();
		double psi = this.getPsi();
		//Log.warning("psi " + psi);
		if (psi > 0) {
			int attemptNr = 0;
	        for (attemptNr = 0; attemptNr < nattempts; attemptNr++) {
	        	
	        	double p = 0;
	        	try {
					p = this.getBlueTreeLogPWithSampling((Tree)tree, lambda, mu, psi);
					if (p == Double.NEGATIVE_INFINITY || p == Double.POSITIVE_INFINITY) throw new Exception();
				} catch (Exception e) {
					double s = Randomizer.nextFloat();
					psi = s*mu / (1-s);
					continue;
				}
	        	
	        	//Log.warning("nr " + attemptNr + " " + lambda + " " + mu + " " + psi + " p=" + p);
	        	
	        	
	        	break;
	        	
	        }
	        
	        if (attemptNr >= nattempts) {
	        	throw new IllegalArgumentException("Cannot find valid initial state. Try tweaking lambda, mu, and psi");
	        }
	        
	        
	        samplingProportionInput.get().setValue(psi / (psi + mu));
	        
	        
		}
        
		
		
		
        
    }
	
	
	@Override
    public double calculateTreeLogLikelihood(final TreeInterface treeInterface) {
        
		
		
		// Get variables
		double p = 0;
		Tree tree = (Tree) treeInterface;
		double lambda = lambdaInput.get().getValue();
		double mu = this.getMu();
		double psi = this.getPsi();
		
		Stubs stubs = stubsInput.get();
		
		
		// Condition checker
		if (lambda <= mu) {
			if (initialising) Log.warning("Cannot initialise because lambda is less than mu");
			initialising = false;
			return Double.NEGATIVE_INFINITY;
		}
		
		
		 // This dist is conditional on ntaxa, so ensure the number does not change
        int n = 0;
		for (Node leaf : tree.getNodesAsArray()) {
			if (leaf.isLeaf() && !leaf.isDirectAncestor() && leaf.getHeight() <= HEIGHT_THRESHOLD_OF_LEAF) {
				n++;
			}
		}
		
		if (nExtantTaxa == -1) {
			nExtantTaxa = n;
		}else if (n != this.nExtantTaxa) {
			return Double.NEGATIVE_INFINITY;
		}
		
		
        
		
		double blueLogP = 0;
		double redLogP = 0;
		
		
		try {
			
			if (!ignoreTreePriorInput.get()) {
				
				// Numerical stability checker
				if (Math.exp(psi) == Double.POSITIVE_INFINITY || Math.exp(-psi) == 0) {
					return Double.NEGATIVE_INFINITY;
				}
				
				if (samplingProportionInput.get() != null) {
					blueLogP = getBlueTreeLogPWithSampling(tree, lambda, mu, psi);
				}
				
				// Psi = 0 - should give the same answer as above when psi=0 but faster
				else {
					blueLogP = getBlueTreeLogPWithoutSampling(tree, lambda, mu);
				}
			
			}
			
			
			
			// Red tree stub density
			if (!ignoreStubPriorInput.get() && stubs != null) {
				
				
				// Do not use the sampling model unless there is at least 1 dated tip
				boolean thereAreDatedTips = false;
				for (Node node : tree.getNodesAsArray()) {
					if (node.isLeaf() && node.getHeight() > 0) {
						thereAreDatedTips = true;
						break;
					}
				}
				
				
				// Confirm that sampled ancestor tips do not have stubs
				for (Node node : tree.getNodesAsArray()) {
					if (node.isDirectAncestor()) {
						int nstubs = stubs.getNStubsOnBranch(node.getNr());
						if (nstubs > 0) {
							//Log.warning("XXXX " + nstubs);
							return Double.NEGATIVE_INFINITY;
						}
					}
				}
				
				
				// Psi > 0
				if (thereAreDatedTips && samplingProportionInput.get() != null && psi > 0) {
					for (int branchNr = 0; branchNr < tree.getNodeCount(); branchNr++) {
						redLogP += getLogPForBranchWithSampling(branchNr, lambda, mu, psi);
					}
					
				}
				
				
				// Psi = 0 - should give the same answer as above when psi=0 but faster
				else {
					for (int branchNr = 0; branchNr < tree.getNodeCount(); branchNr++) {
						redLogP += getLogPForBranch(branchNr, lambda, mu);
					}
					
				}
				
			}
			
		
		}catch (Exception e) {
			
			//Log.warning("numerical");
			//e.printStackTrace();
			
			// Numerical errors
			return Double.NEGATIVE_INFINITY;
		}
		
		
		if (initialising && blueLogP == Double.NEGATIVE_INFINITY) {
			Log.warning("Cannot initialise stumped tree prior because blue p is negative infinity");
		}
		
		else if (initialising && redLogP == Double.NEGATIVE_INFINITY) {
			Log.warning("Cannot initialise stumped tree prior because red p is negative infinity");
		}
		
		
		if (p == Double.POSITIVE_INFINITY) {
			return Double.NEGATIVE_INFINITY;
		}
		
		
		p = blueLogP + redLogP;
		initialising = false;
		return p;
		
		
		
    }
	
	
	
	// Assumes that psi=0, rho=1
	public double getBlueTreeLogPWithoutSampling(Tree tree, double lambda, double mu) {
		
		
		// Blue tree speciation density, conditional on survival of both children
		double blueLogP = -this.getBlueIntegral();
		for (int branchNr = 0; branchNr < tree.getNodeCount(); branchNr++) {
			if (tree.getNode(branchNr).isLeaf()) continue;
			double height = tree.getNode(branchNr).getHeight();
			double logRate = this.getBlueTreeLogBirthRate(height, lambda, mu);
			blueLogP += logRate;
		}
		
		return blueLogP;
		
	}
	
	
	// Equation 9 from
	// Stadler, Tanja. "Sampling-through-time in birth–death trees." Journal of theoretical biology 267.3 (2010): 396-404.
	public double getBlueTreeLogPWithSampling(Tree tree, double lambda, double mu, double psi) throws Exception {
		
		
		
		double rho=1;
		
		double m = 0, n = 0, k = 0;
		for (Node leaf : tree.getNodesAsArray()) {
			
			if (!leaf.isLeaf()) continue;
			
			
			// Sampled ancestor
			if (leaf.isDirectAncestor()) {
				k++;
			}
			
			// Extant leaf
			else if (!leaf.isDirectAncestor() && leaf.getHeight() <= HEIGHT_THRESHOLD_OF_LEAF) {
				n++;
			}
			
			// Extinct leaf
			else if (!leaf.isDirectAncestor() && leaf.getLength() > HEIGHT_THRESHOLD_OF_LEAF) {
				m++;
			}
			
			else {
				//Log.warning("Dev Error 3543: " + leaf.getID() +" is not a valid leaf " + leaf.getLength() + " " + leaf.isDirectAncestor());
			}
			
			
			
		}
		
		if (m + n + k != tree.getLeafNodeCount()) {
			//Log.warning("n=" + n + " m=" + m + " k=" + k + " != " + tree.getLeafNodeCount());
			return Double.NEGATIVE_INFINITY;
		}
		
		
		double blueLogP = 0;
		
		// Equation 9: p(T|n)
		double c1 = getC1(lambda, mu, psi, rho); 
		double c2 = getC2(lambda, mu, psi, rho); 
		double x1 = tree.getRoot().getHeight();
		blueLogP += Math.log(4) + Math.log(n) + Math.log(rho) + Math.log(lambda) + Math.log(psi)*(k+m); // maybe extra lambda is typo?? 
		blueLogP -= Math.log(c1) + Math.log(c2+1) + Math.log(1 - c2 + (1+c2)*Math.exp(c1*x1)); // x=x1, typo in paper
		
		
		// Loop through internal nodes including root
		for (Node internal : tree.getNodesAsArray()) {
			
			if (internal.isLeaf()) continue;
			
			// Confirm that neither child is a sampled ancestor
			if (internal.isFake()) continue;
			
			
			double x = internal.getHeight();
			blueLogP += Math.log(lambda) + Math.log(getP1(x, lambda, mu, psi, rho, c1, c2));
		}
		
		// Loop through extinct leaves
		for (Node leaf : tree.getNodesAsArray()) {
			if (!leaf.isLeaf()) continue;
			if (leaf.getHeight() <= HEIGHT_THRESHOLD_OF_LEAF) continue;
			if (leaf.isDirectAncestor()) continue;
			double y = leaf.getHeight();
			blueLogP += Math.log(getP0(y, lambda, mu, psi, rho, c1, c2)) - Math.log(getP1(y, lambda, mu, psi, rho, c1, c2));
		}
		

		
		return blueLogP;
	}
	
	
	private double getC1(double lambda, double mu, double psi, double rho) {
		return Math.sqrt(Math.abs(Math.pow(lambda - mu - psi, 2) + 4*lambda*psi));
	}
	
	private double getC2(double lambda, double mu, double psi, double rho) {
		double c1 = getC1(lambda, mu, psi, rho); 
		return -(lambda-mu-2*lambda*rho-psi) / c1;
	}
	
	// Equation 1 from
	// Stadler, Tanja. "Sampling-through-time in birth–death trees." Journal of theoretical biology 267.3 (2010): 396-404.
	private double getP0(double height, double lambda, double mu, double psi, double rho, double c1, double c2) throws Exception {
		double a = lambda + mu + psi;
		double exp = Math.exp(-c1*height)*(1-c2);
		double b = 1+c2;
		double top = a + c1*(exp - b)/(exp + b);
		double bottom = 2*lambda;
		double result = top/bottom;
		if (Double.isNaN(result) || result <= 0) {
			//Log.warning("c1=" + c1 + " c2=" + c2 + " lambda =" + lambda + " mu=" + mu + " psi=" + psi);
			throw new Exception("Numerical error: p0 is non-positive " + top + "/" + bottom + "=" + result);
		}
		return result;
	}
	
	// Equation 2 from
	// Stadler, Tanja. "Sampling-through-time in birth–death trees." Journal of theoretical biology 267.3 (2010): 396-404.
	private double getP1(double height, double lambda, double mu, double psi, double rho, double c1, double c2) throws Exception {
		double top = 4*rho;
		double bottom = 2*(1-c2*c2) + Math.exp(-c1*height)*(1-c2)*(1-c2) + Math.exp(c1*height)*(1+c2)*(1+c2);
		double result = top/bottom;
		if (Double.isNaN(result) || result <= 0) {
			//Log.warning("c1=" + c1 + " c2=" + c2 + " lambda =" + lambda + " mu=" + mu + " psi=" + psi);
			throw new Exception("Numerical error: p1 is non-positive " + top + "/" + bottom + "=" + result);
		}
		return result;
	}
	
	
	// Theorem 3.3 from
	// Stadler, Tanja. "Sampling-through-time in birth–death trees." Journal of theoretical biology 267.3 (2010): 396-404.
	private double getP0Hat(double t, double lambda, double mu, double rho) throws Exception {
		
		double e = Math.exp(-(lambda-mu)*t);
		double a = rho*(lambda - mu);
		double b = rho*lambda + (lambda*(1-rho) - mu)*e;
		return 1 - (a/b);
		
	}
	
	
	// Theorem 3.3 from
	// Stadler, Tanja. "Sampling-through-time in birth–death trees." Journal of theoretical biology 267.3 (2010): 396-404.
	private double getP1Hat(double t, double lambda, double mu, double rho) throws Exception {
		
		double e = Math.exp(-(lambda-mu)*t);
		double a = rho*(lambda - mu)*(lambda - mu)*e;
		double b = rho*lambda + (lambda*(1-rho) - mu)*e;
		return a/(b*b);
		
	}
	
	
	
	/**
	 * Red tree with psi>0
	 * @param branchNr
	 * @param lambda
	 * @param mu
	 * @return
	 */
	private double getLogPForBranchWithSampling(int branchNr, double lambda, double mu, double psi) {
		
		
		
		double p = 0;
		double rho = 1;
		
		Tree tree = (Tree) treeInput.get();
		Stubs stubs = stubsInput.get();
		double T = tree.getRoot().getHeight();
		Node node = tree.getNode(branchNr);
		
		if (!stubs.estimateStubs()) return 0;

		
		
		// Reverse times
		double h0 = node.isRoot() ? T : node.getParent().getHeight();
		double h1 = node.getHeight();
		
		
		// Just count the stubs - don't care about time
		if (!stubs.getReversibleJump()) {
			
			
			
			if (node.isRoot()) return 0;
			
			// A zero length branch must have zero stubs
			int m = stubs.getNStubsOnBranch(branchNr);
			if (h0 <= h1) {
				if (m == 0) return 0;
				return Double.NEGATIVE_INFINITY;
			}
			
			
			double g = getGForInterval(h1, h0, lambda, mu, psi, rho);
			
			
			// Poisson(g) distribution
			p = m*Math.log(g) - g;
			for (int i = 2; i <= m; i ++) p += -Math.log(i);
			
			//Log.warning("g=" + g);
			//Log.warning( " m=" + m + " " + g + " " + p + " " + h0 + " " + h1);
			
			return p;
			
		}
		
		// Stubs on this branch
		int m = 0;
		for (int stubNr = 0; stubNr < stubs.getStubDimension(); stubNr++) {
			if (!stubs.includeStub(stubNr)) continue;
			if (stubs.getBranches().getValue(stubNr) == branchNr) {
				double stubHeight = stubs.getAbsoluteTimeOfStub(stubNr);
				
				if (stubHeight < node.getHeight() || stubHeight > node.getParent().getHeight()) {
					return Double.NEGATIVE_INFINITY;
				}
				if (node.isDirectAncestor()) {
					return Double.NEGATIVE_INFINITY;
				}
				
				double logRate = this.getRedTreeBirthLogRate(stubHeight, lambda, mu, psi, rho);
				p += logRate;
				m++;
			}
		}
		
		// Integral and combinatorial term across the branch
		double g = getGForInterval(h1, h0, lambda, mu, psi, rho);
		//Log.warning("g=" + g);
		p += -g;
		for (int i = 2; i <= m; i ++) p += -Math.log(i);
		
		return p;
	}
	
	
	/**
	 * Red tree 
	 * @param branchNr
	 * @return
	 */
	private double getLogPForBranch(int branchNr, double lambda, double mu) {
		

		// Get variables
		double p = 0;
		Tree tree = (Tree) treeInput.get();
		Stubs stubs = stubsInput.get();
		double T = tree.getRoot().getHeight();
		
		if (!stubs.estimateStubs()) return 0;
		
		
		Node node = tree.getNode(branchNr);
		
		// Reverse times
		double h0 = node.isRoot() ? T : node.getParent().getHeight();
		double h1 = node.getHeight();
		
		
		// Just count the stubs - don't care about time
		if (!stubs.getReversibleJump()) {
			
			
			if (node.isRoot()) return 0;
			
			// A zero length branch must have zero stubs
			int m = stubs.getNStubsOnBranch(branchNr);
			if (h0 <= h1) {
				if (m == 0) return 0;
				return Double.NEGATIVE_INFINITY;
			}
			
			double g = getGForInterval(h1, h0, lambda, mu);
			
			// Poisson(g) distribution
			p = m*Math.log(g) - g;
			for (int i = 2; i <= m; i ++) p += -Math.log(i);
			return p;
			
		}
		
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
		double g = getGForInterval(h1, h0, lambda, mu);
		p += -g;
		for (int i = 2; i <= m; i ++) p += -Math.log(i);
		
		return p;
		
	}
	
	public double getPsi() {
		if (samplingProportionInput.get() == null) return 0.0;
		double mu = this.getMu();
		double samplingProportion = samplingProportionInput.get().getValue();
		return samplingProportion*mu / (1-samplingProportion);
	}
	

	
	
	public double getMu() {
		double lambda = lambdaInput.get().getValue();
		if (r0Input.get() != null) {
			return lambda / r0Input.get().getValue();
		}else {
			return turnoverInput.get().getValue() * lambda;
		}
	}
	
	

	
	
	private double getGForInterval(double h1, double h0, double lambda, double mu) {
		
		Tree tree = (Tree) treeInput.get();
		double T = tree.getRoot().getHeight();
		
		return getRedRateIntegral(h1, lambda, mu, T) - getRedRateIntegral(h0, lambda, mu, T);
	}
	
	// Equation 1 from
	// https://www.biorxiv.org/content/10.1101/2024.04.15.589468v1.full.pdf
	private double getGForInterval(double h1, double h0, double lambda, double mu, double psi, double rho) {
		
		Tree tree = (Tree) treeInput.get();
		double T = tree.getRoot().getHeight();
		
		
		double c1 = getC1(lambda, mu, psi, rho); 
		double c2 = getC2(lambda, mu, psi, rho); 
		double f = ((c2-1)*Math.exp(-c1*h1) -c2 - 1) / ((c2-1)*Math.exp(-c1*h0) -c2 - 1);
		return (h0-h1) * (lambda + mu + psi -c1) + 2*Math.log(f); 
		
		
		//return getRedRateIntegral(h1, lambda, mu, psi, T) - getRedRateIntegral(h0, lambda, mu, psi, T);
	}

	
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

	

	

	
	
	private double getGSum() {
		
		double gSum = 0;
		
		Tree tree = (Tree) treeInput.get();
		double lambda = lambdaInput.get().getValue();
		double mu = getMu();
		double T = tree.getRoot().getHeight();
		
		double psi = this.getPsi();
		double rho = 1;
		
		// For each branch
		for (int branchNr = 0; branchNr < tree.getNodeCount(); branchNr++) {
			
			
			Node node = tree.getNode(branchNr);
			
			// Reverse times
			double h0 = node.isRoot() ? (node.getHeight()) : node.getParent().getHeight();
			double h1 = node.getHeight();
			
			// Integral term across the branch
			double g;
			if (!ignoreStubPriorInput.get()) {
				g = getGForInterval(h1, h0, lambda, mu, psi, rho);
			}else {
				g = getRedRateIntegral(h1, lambda, mu, T) - getRedRateIntegral(h0, lambda, mu, T);
			}
			
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
		double mu = getMu();
		double T = tree.getRoot().getHeight();
		
		
		double g = 0;
		for (Node node : tree.getNodesAsArray()) {
			
			
			// Reverse times
			double h0 = node.isRoot() ? (node.getHeight()) : node.getParent().getHeight();
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
		double mu = getMu();
		double T = tree.getRoot().getHeight();
		
		
		double blueIntegral = 0;
		for (Node node : tree.getNodesAsArray()) {
			
			
			// Reverse times
			double h0 = node.isRoot() ? (node.getHeight()) : node.getParent().getHeight();
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
	
	
	private double getRedTreeBirthLogRate(double timeRemaining, double lambda, double mu, double psi, double rho) {
		double logqt = this.getLogQt(timeRemaining, lambda, mu, psi, rho);
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
	
	
	private double getLogQt(double time, double lambda, double mu, double psi, double rho) {
		
		double c1 = getC1(lambda, mu, psi, rho); 
		double c2 = getC2(lambda, mu, psi, rho); 
		
		
		double exp = Math.exp(-c1*time);
		double top = exp*(1-c2) - (1+c2);
		double btm = exp*(1-c2) + (1+c2);
		double p = lambda + mu + psi + c1 * top/btm;
		return Math.log(p) - Math.log(2*lambda);
	}
	
	

	
	@Override
    protected boolean requiresRecalculation() {
		
		
		if (pInput.get() != null && pInput.get().somethingIsDirty()) {
			return true;
		}
		
		
        return super.requiresRecalculation() || 
        		InputUtil.isDirty(lambdaInput) || 
        		InputUtil.isDirty(r0Input) || 
        		//InputUtil.isDirty(rhoInput) || 
        		InputUtil.isDirty(turnoverInput) || 
        		InputUtil.isDirty(samplingProportionInput) ||
        		InputUtil.isDirty(stubsInput);
    }

	
    /**
     * Loggable interface implementation follows *
     */
    @Override
    public void init(final PrintStream out) {
       // out.print(getID() + "\t" + "g0\tg1\tg2\tg3\t"+ "stubGsum\t"+ "stubH\t");
        out.print(getID() + "\t");
    }

    @Override
    public void log(final long sample, final PrintStream out) {
    	//out.print(getCurrentLogP() + "\t" + this.getG(0) + "\t" + this.getG(1) + "\t" + this.getG(2) + "\t" + this.getG(3) + "\t" + this.getGSum()  +"\t"+ this.getH() + "\t");
        out.print(getCurrentLogP() + "\t");
    }


	public double getTotalTime() {
		Tree tree = (Tree) treeInput.get();
		return tree.getRoot().getHeight();
	}


	public Tree getTree() {
		return (Tree) treeInput.get();
	}


	@Override
	public double getMeanStubNumber(double height, double parentHeight) {
		
		if (parentHeight <= height) return 0;
		
		double lambda = lambdaInput.get().getValue();
		double mu = this.getMu();
		double psi = this.getPsi();
		if (psi == 0) {
			return getGForInterval(height, parentHeight, lambda, mu);
		}else {
			return getGForInterval(height, parentHeight, lambda, mu, psi, 1);
		}
	}

	
	



	
	
}


