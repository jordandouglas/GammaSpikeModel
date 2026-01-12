package gammaspike.distribution;


import java.io.PrintStream;

import beast.base.core.*;
import beast.base.core.Input.Validate;
import beast.base.evolution.speciation.SpeciesTreeDistribution;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import gammaspike.tree.Stubs;



@Description("Prior distribution on a stumped tree")
@Citation(value =
"Douglas, J., Bouckaert, R., Harris, S.C., Carter Jr, C.W., Wills, P.R. (2025) Evolution is coupled with branching across many granularities of life. Proceedings of the Royal Society Series B 29220250182", DOI = "http://doi.org/10.1098/rspb.2025.0182",
year = 2025, firstAuthorSurname = "Douglas")
public class StumpedTreePrior extends SpeciesTreeDistribution implements StubExpectation {

	final public Input<RealParameter> originInput = new Input<>("origin", "the time when the process started", Validate.OPTIONAL);
	// if the tree likelihood is condition on sampling at least one individual then set to true one of the inputs:
	public Input<Boolean> conditionOnSamplingInput = new Input<Boolean>("conditionOnSampling", "the tree " +
			"likelihood is conditioned on sampling at least one individual if condition on origin", false);
	public Input<Boolean> conditionOnRhoSamplingInput = new Input<Boolean>("conditionOnRhoSampling", "the tree " +
			"likelihood is conditioned on sampling at least one individual in present if condition on origin", false);

	public Input<Boolean> conditionOnRootInput = new Input<Boolean>("conditionOnRoot", "the tree " +
			"likelihood is conditioned on the root height otherwise on the time of origin", false);

	final public Input<Function> lambdaInput = new Input<>("lambda", "Birth rate (lambda)", Validate.REQUIRED);
	final public Input<Function> r0Input = new Input<>("r0", "reproduction number, or lambda / mu", Validate.REQUIRED);
	final public Input<Function> psiInput = new Input<>("psi", "Sampling rate (psi)", Validate.OPTIONAL);

	final public Input<Function> netDiversificationRateInput = new Input<>("netDiversificationRate", "net diversification rate, lambda - mu", Validate.XOR, lambdaInput);
	final public Input<Function> turnoverInput = new Input<>("turnover", "the ratio of death (mu) and birth (lambda) rates, mu / lambda", Validate.XOR, r0Input);
	final public Input<Function> samplingProportionInput = new Input<>("samplingProportion", "sampling rate psi divided by (psi + mu)", Validate.OPTIONAL);

//	final public Input<Function> r0Input = new Input<>("r0", "the ratio of birth (lambda) and death (mu) rates, lambda / mu");

	final public Input<RealParameter> rhoInput = new Input<>("rho", "extant sampling probability (rho)", Validate.OPTIONAL);

	final public Input<Stubs> stubsInput = new Input<>("stubs", "the stubs of this tree, if the stub-free inference approach (integrating over the number of stubs) is used, this is ignored", Input.Validate.OPTIONAL);


	final public Input<Boolean> ignoreTreePriorInput = new Input<>("ignoreTreePrior", "ignore tree prior (debugging)", false);
	final public Input<Boolean> ignoreStubPriorInput = new Input<>("ignoreStubPrior", "ignore tree prior (debugging)", false);

//	final public Input<RealParameter> pInput = new Input<>("p", "the probability of a single-mutant. if included, then stubs correspond to bursts, and"
//			+ " the rate of stub appearance is multiplied by 1-p/2", Input.Validate.OPTIONAL);
//	final public Input<Boolean> perBranchSpikeInput = new Input<>("perBranchSpike", "if true, then the indicator will disregard labels", false);


	final double HEIGHT_THRESHOLD_OF_LEAF = 1e-10;
	int nExtantTaxa; // n

	boolean initialising = true;
	protected boolean originSpecified;


	@Override
	public boolean canHandleTipDates() {
		return samplingProportionInput.get() != null;
	}

	@Override
    public void initAndValidate() {
        super.initAndValidate();

		originSpecified = originInput.get() != null;

//		if (originSpecified && !conditionOnSamplingInput.get() && !conditionOnRhoSamplingInput.get()) {
//			throw new IllegalArgumentException("Either conditionOnSampling or conditionOnRhoSampling has to be set to \"true\" when an origin time is specified.");
//		}
//		if (!originSpecified && (conditionOnSamplingInput.get() || conditionOnRhoSamplingInput.get())) {
//			throw new IllegalArgumentException("Cannot use conditionOnSampling or conditionOnRhoSampling without specifying an origin time.");
//		}
//		if (!originSpecified && !conditionOnRootInput.get()) {
//			throw new IllegalArgumentException("Specify origin or set conditionOnRoot input to \"true\"");
//		}
		
		if (psiInput.get() != null && samplingProportionInput.get() != null) {
			throw new RuntimeException("Please specifiy either psi or samplingProportion or neither, but not both");
		}
		
		if (originSpecified && conditionOnRootInput.get()) {
			throw new RuntimeException("Remove origin or set conditionOnRoot input to \"false\"");
		}
		if (conditionOnSamplingInput.get() && conditionOnRhoSamplingInput.get()){
			throw new IllegalArgumentException("Either set to \"true\" only one of conditionOnSampling and conditionOnRhoSampling inputs or don't specify both!");
		}
		if (conditionOnRootInput.get() && !conditionOnRhoSamplingInput.get() && !conditionOnSamplingInput.get()) {
			throw new IllegalArgumentException("When conditioning on the root, we always assume that both sides of the initial bifurcation event are sampled. Please set either " +
					"conditionOnSampling or conditionOnRhoSampling to true.");
		}

       TreeInterface tree = treeInput.get();
        if (tree == null) {
            tree = treeIntervalsInput.get().treeInput.get();
        }

		double rootHeight = tree.getRoot().getHeight();
		if (originSpecified && origin() < rootHeight){
			throw new IllegalArgumentException("Initial value of origin (" + origin() + ") should be greater than initial root height (" +rootHeight + ")");
		}

		// This distribution is conditional on the number of extant taxa n, so ensure the number does not change
        nExtantTaxa = -1;

        // Ensure valid initial state
//        final int MAX_ATTEMPTS = 1000;
//		double lambda = this.getLambda();
//		double mu = this.getMu();
//		double psi = this.getPsi();
//		double rho = this.getRho();
		//Log.warning("psi " + psi);
//		if (psi > 0 && (samplingProportionInput.get() instanceof RealParameter)) {
//            int attemptNr;
//            for (attemptNr = 0; attemptNr < MAX_ATTEMPTS; attemptNr++) {
//
//                double p = 0;
//                try {
//                    p = this.getBlueTreeLogPWithSampling((Tree) tree, lambda, mu, psi, rho);
//                    if (p == Double.NEGATIVE_INFINITY || p == Double.POSITIVE_INFINITY) throw new Exception();
//                } catch (Exception e) {
//                    double s = Randomizer.nextFloat();
//                    psi = s * mu / (1 - s);
//                    continue;
//                }
//
//                Log.warning("nr " + attemptNr + " " + lambda + " " + mu + " " + psi + " p=" + p);
//                break;
//
//            }
//
//            if (attemptNr >= MAX_ATTEMPTS) {
//                throw new IllegalArgumentException("Cannot find valid initial state. Try tweaking lambda, mu, psi, and rho");
//            }
//
//            ((RealParameter) samplingProportionInput.get()).setValue(psi / (psi + mu));
//        }

    }


	@Override
    public double calculateTreeLogLikelihood(final TreeInterface treeInterface) {
		Tree tree = (Tree) treeInterface;
		double lambda = this.getLambda();
		double mu = this.getMu();
		double psi = this.getPsi();
		double rho = this.getRho();

		Stubs stubs = stubsInput.get();

		// Condition checker
		if (lambda <= mu) {
			if (initialising) Log.warning("Cannot initialise because lambda is less than mu");
			initialising = false;
			return Double.NEGATIVE_INFINITY;
		}

		double logP = 0;
		double treeLogP = 0; // Probability of observing tree T conditional on number of extant taxa n
		double stubLogP = 0; // Probability of observing B stubs at ages Z

		try {

			// Tree density
			// Probability of observing tree T conditional on number of extant taxa n
			if (!ignoreTreePriorInput.get()) {
				// Numerical stability checker
				if (Math.exp(psi) == Double.POSITIVE_INFINITY || Math.exp(-psi) == 0) {
					return Double.NEGATIVE_INFINITY;
				}
				if (originSpecified || conditionOnRootInput.get()) {
					double x0 = originSpecified ? origin() : Double.POSITIVE_INFINITY;
					double x1 = tree.getRoot().getHeight();
					if (x0 < x1) {
						return Double.NEGATIVE_INFINITY;
					}
					treeLogP = getTreeLogPConditionOnOriginOrRoot(tree, lambda, mu, psi, rho, x0, x1);
				}
				else {

					// This distribution is conditional on the number of extant taxa n, so ensure the number does not change
					int n = 0;
					for (Node leaf : tree.getNodesAsArray()) {
						if (leaf.isLeaf() && !leaf.isDirectAncestor() && leaf.getHeight() <= HEIGHT_THRESHOLD_OF_LEAF) {
							n++;
						}
					}
					if (nExtantTaxa == -1) {
						nExtantTaxa = n;
					} else if (n != this.nExtantTaxa) {
						return Double.NEGATIVE_INFINITY;
					}

					// Tree likelihood calculation when psi = 0 and rho = 1
					if ((psi == 0) && (rho == 1)) {
						treeLogP = getTreeLogPWithoutPsiSamplingCompleteRho(tree, lambda, mu);
					}
					// Tree likelihood calculation when psi = 0 and rho < 1
					else if (psi == 0) {
						treeLogP = getTreeLogPWithoutPsiSampling(tree, lambda, mu, rho);
					}
					// Tree likelihood calculation when psi > 0 and rho < 1
					else {
						treeLogP = getTreeLogPWithPsiSampling(tree, lambda, mu, psi, rho);
					}
				}
			}

			// Stub density
			// Probability of observing B stubs at ages Z
			if (!ignoreStubPriorInput.get() && stubs != null) {

				// Do not use the sampling model unless there is at least 1 dated tip
				// A dated tip means non-extant sampling
//				boolean thereAreDatedTips = false;
//				for (Node node : tree.getNodesAsArray()) {
//					if (node.isLeaf() && node.getHeight() > 0) {
//						thereAreDatedTips = true;
//						break;
//					}
//				}

				// Confirm that sampled ancestor tips do not have stubs
				for (Node node : tree.getNodesAsArray()) {
					if (node.isDirectAncestor()) {
						int stubCount = stubs.getNStubsOnBranch(node.getNr());
						if (stubCount > 0) {
							//Log.warning("XXXX " + nstubs);
							return Double.NEGATIVE_INFINITY;
						}
					}
				}

				// Stub probability calculation when psi = 0 and rho = 1
				if ((psi == 0) && (rho == 1)) {
					for (int branchNr = 0; branchNr < tree.getNodeCount(); branchNr++) {
						stubLogP += getStubLogPForBranchWithoutSampling(branchNr, lambda, mu);
					}
				}
				// Stub probability calculation when psi > 0 and rho < 1
				else {
					for (int branchNr = 0; branchNr < tree.getNodeCount(); branchNr++) {
						stubLogP += getStubLogPForBranchWithSampling(branchNr, lambda, mu, psi, rho);
					}
				}

			}


		} catch (Exception e) {

			//Log.warning("numerical");
			//e.printStackTrace();

			// Numerical errors
			return Double.NEGATIVE_INFINITY;
		}


		if (initialising && treeLogP == Double.NEGATIVE_INFINITY) {
			Log.warning("Cannot initialise stumped tree prior because tree log p is negative infinity");
		}

		else if (initialising && stubLogP == Double.NEGATIVE_INFINITY) {
			Log.warning("Cannot initialise stumped tree prior because stub log p is negative infinity");
		}

		if (logP == Double.POSITIVE_INFINITY) {
			return Double.NEGATIVE_INFINITY;
		}

		logP = treeLogP + stubLogP;
		initialising = false;
		return logP;

    }


	// Equation 9 from
	// Stadler, Tanja. "Sampling-through-time in birth–death trees." Journal of theoretical biology 267.3 (2010): 396-404.
	// Tree likelihood calculation when psi > 0 and rho < 1
	public double getTreeLogPWithPsiSampling(Tree tree, double lambda, double mu, double psi, double rho) throws Exception {
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
			//else {
				//Log.warning("Dev Error 3543: " + leaf.getID() +" is not a valid leaf " + leaf.getLength() + " " + leaf.isDirectAncestor());
			//}

		}
		if (m + n + k != tree.getLeafNodeCount()) {
			//Log.warning("n=" + n + " m=" + m + " k=" + k + " != " + tree.getLeafNodeCount());
			return Double.NEGATIVE_INFINITY;
		}

		double treeLogP = 0;

		// Equation 9: p(T|n)
		// Conditioned on sampling n extant individuals
		double c1 = getC1(lambda, mu, psi);
		double c2 = getC2(lambda, mu, psi, rho);
		double x1 = tree.getRoot().getHeight();
		treeLogP += Math.log(4) + Math.log(n) + Math.log(rho) + Math.log(lambda) + Math.log(psi) * (k + m);
		treeLogP -= Math.log(c1) + Math.log(c2+1) + Math.log(1 - c2 + (1 + c2) * Math.exp(c1 * x1));

		// Loop through internal nodes including root
		for (Node internal : tree.getNodesAsArray()) {
			if (internal.isLeaf()) continue;
			if (internal.isFake()) continue; // Confirm that neither child is a sampled ancestor
			double nodeHeight = internal.getHeight();
			treeLogP += Math.log(lambda) + Math.log(getP1(nodeHeight, rho, c1, c2));
		}

		// Loop through extinct leaves
		for (Node leaf : tree.getNodesAsArray()) {
			if (!leaf.isLeaf()) continue;
			if (leaf.getHeight() <= HEIGHT_THRESHOLD_OF_LEAF) continue;
			if (leaf.isDirectAncestor()) continue;
			double leafHeight = leaf.getHeight();
			treeLogP += Math.log(getP0(leafHeight, lambda, mu, psi, c1, c2)) - Math.log(getP1(leafHeight, rho, c1, c2));
		}

		return treeLogP;

	}


	// Tree likelihood calculation when psi = 0 and rho < 1
	public double getTreeLogPWithoutPsiSampling(Tree tree, double lambda, double mu, double rho) throws Exception {
		double n = 0;
		for (Node leaf : tree.getNodesAsArray()) {
			if (!leaf.isLeaf()) continue;
			// Extant leaf
			if (!leaf.isDirectAncestor() && leaf.getHeight() <= HEIGHT_THRESHOLD_OF_LEAF) {
				n++;
			}
		}
		if (n != tree.getLeafNodeCount()) {
			return Double.NEGATIVE_INFINITY;
		}

		double treeLogP = 0;

		double x1 = tree.getRoot().getHeight();
		treeLogP += Math.log(n) + Math.log(lambda - mu) - (lambda - mu) * x1;
		treeLogP -= Math.log(rho * lambda + (lambda * (1 - rho) - mu) * Math.exp(-(lambda - mu) * x1));

		// Loop through internal nodes including root
		for (Node internal : tree.getNodesAsArray()) {
			if (internal.isLeaf()) continue;
			if (internal.isFake()) continue; // Confirm that neither child is a sampled ancestor
			double nodeHeight = internal.getHeight();
			treeLogP += Math.log(lambda) + Math.log(rho) + 2 * Math.log(lambda - mu) - (lambda - mu) * nodeHeight;
			treeLogP -= 2 * Math.log(rho * lambda + (lambda * (1 - rho) - mu) * Math.exp(-(lambda - mu) * nodeHeight));
		}
		return treeLogP;
	}


	// Tree likelihood calculation when psi = 0 and rho = 1
	public double getTreeLogPWithoutPsiSamplingCompleteRho(Tree tree, double lambda, double mu) {
		// Tree probability density, conditional on survival of both children
		double treeLogP = -this.getBlueIntegral(tree, lambda, mu); // Log-likelihood of the tree surviving from its root to the present day without any speciation events occurring within the internal branches
		// Loop through all branches (nodes)
		for (int branchNr = 0; branchNr < tree.getNodeCount(); branchNr++) {
			// Skip leaf nodes (tips)
			if (tree.getNode(branchNr).isLeaf()) continue;
			double nodeHeight = tree.getNode(branchNr).getHeight();
			double logRate = this.getTreeLogBirthRate(nodeHeight, lambda, mu); // Sum of log-likelihood contribution of each speciation event (internal node)
			treeLogP += logRate;
		}
		return treeLogP;
	}

	// Methods used by tree likelihood calculation when psi = 0 and rho = 1
	private double getBlueIntegral(Tree tree, double lambda, double mu) {
		double treeHeight = tree.getRoot().getHeight();

		double blueIntegral = 0;
		// Loop through all nodes (all branches)
		for (Node node : tree.getNodesAsArray()) {

			// Get the origin time h0 and end time h1 of each branch branchNr
			double h0 = node.isRoot() ? (node.getHeight()) : node.getParent().getHeight();
			double h1 = node.getHeight();
			// Sum of log-likelihood of each branch's survival without splitting
			double integralBranch = getBlueRateIntegral(h1, lambda, mu, treeHeight) - getBlueRateIntegral(h0, lambda, mu, treeHeight); // Integral over a single branch interval between h1 and h0
			blueIntegral += integralBranch;

		}

		return blueIntegral;

	}
	private double getBlueRateIntegral(double height, double lambda, double mu, double rootHeight) {
		double t = rootHeight - height; // Time elapsed since the root
		double tlm = t * (lambda - mu);
		double b = Math.log(lambda * Math.exp(lambda * rootHeight) - mu * Math.exp(tlm + mu * rootHeight));
		return tlm - b;
	}
	private double getTreeLogBirthRate(double timeRemaining, double lambda, double mu) {
		double logQt = this.getLogQt(timeRemaining, lambda, mu);
		return Math.log(lambda) + Math.log(1 - Math.exp(logQt)); // lambda * (1 - p0(t)) birth rate conditioned on survival to the present day
	}

	public double getTreeLogPConditionOnOriginOrRoot(Tree tree, double lambda, double mu, double psi, double rho, double x0, double x1) {
		double treeLogP;
		int nodeCount = tree.getNodeCount();
		double c1 = getC1(lambda, mu, psi);
		double c2 = getC2(lambda, mu, psi, rho);

		if (!conditionOnRootInput.get()){
			treeLogP = -Math.log(q(x0, c1, c2));
		} else {
			if (tree.getRoot().isFake()){ //when conditioning on the root we assume the process
				//starts at the time of the first branching event and
				//that means that the root can not be a sampled ancestor
				return Double.NEGATIVE_INFINITY;
			} else {
				treeLogP = -Math.log(q(x1, c1, c2));
			}
		}

		if (conditionOnSamplingInput.get()) {
			if (conditionOnRootInput.get()) {
				treeLogP -= Math.log(lambda * (1 - getP0(x1, lambda, mu, psi, c1, c2)) * (1 - getP0(x1, lambda, mu, psi, c1, c2)));
			} else {
				treeLogP -= Math.log(1 - getP0(x0, lambda, mu, psi, c1, c2));
			}

		}

		if (conditionOnRhoSamplingInput.get()) {
			if (conditionOnRootInput.get()) {
				treeLogP -= Math.log(lambda*oneMinusP0Hat(x1, lambda, mu, rho)* oneMinusP0Hat(x1, lambda, mu, rho));
			}  else {
				treeLogP -= Math.log(oneMinusP0Hat(x0, lambda, mu, rho));
			}
		}

		int internalNodeCount = tree.getLeafNodeCount() - tree.getDirectAncestorNodeCount() - 1;

		treeLogP += internalNodeCount*Math.log(2);

		for (int i = 0; i < nodeCount; i++) {
			if (tree.getNode(i).isLeaf()) {
				if  (!tree.getNode(i).isDirectAncestor())  {
					if (tree.getNode(i).getHeight() > 0.000000000005 || rho == 0.) {
						treeLogP += Math.log(psi) + Math.log(q(tree.getNode(i).getHeight(), c1, c2)) + Math.log(getP0(tree.getNode(i).getHeight(), lambda, mu, psi, c1, c2));
					} else {
						treeLogP += Math.log(4*rho);
					}
				}
			} else {
				if (tree.getNode(i).isFake()) {
					treeLogP += Math.log(psi); // removal probability assumed to be zero, log(1 - r) = 0
				} else {
					treeLogP += Math.log(lambda) - Math.log(q(tree.getNode(i).getHeight(), c1, c2));
				}
			}
		}

		return treeLogP;

	}

	protected double q(double t, double c1, double c2) {
		return Math.exp(c1 * t) * (1 + c2) * (1 + c2) + Math.exp(-c1 * t) * (1 - c2) * (1 - c2) + 2 * (1 - c2 * c2);
	}
	protected double oneMinusP0Hat(double t, double lambda, double mu, double rho) {
		return rho*(lambda-mu)/(lambda*rho + (lambda*(1-rho) - mu)* Math.exp((mu-lambda) * t)) ;
	}

	private double getC1(double lambda, double mu, double psi) {
		return Math.sqrt(Math.abs(Math.pow(lambda - mu - psi, 2) + 4 * lambda * psi));
	}

	private double getC2(double lambda, double mu, double psi, double rho) {
		return -(lambda - mu - 2 * lambda * rho - psi) / getC1(lambda, mu, psi);
	}

	// Equation 1 from
	// Stadler, Tanja. "Sampling-through-time in birth–death trees." Journal of theoretical biology 267.3 (2010): 396-404.
	private double getP0(double height, double lambda, double mu, double psi, double c1, double c2) {
		return (lambda + mu + psi - c1 * ((1 + c2) - Math.exp(-c1 * height) * (1 - c2)) / ((1 + c2) + Math.exp(-c1 * height) * (1 - c2))) / (2 * lambda);
//		double a = lambda + mu + psi;
//		double exp = Math.exp(-c1*height)*(1-c2);
//		double b = 1+c2;
//		double top = a + c1*(exp - b)/(exp + b);
//		double bottom = 2*lambda;
//		double result = top/bottom;
//		if (Double.isNaN(result) || result <= 0) {
//			Log.warning("c1=" + c1 + " c2=" + c2 + " lambda =" + lambda + " mu=" + mu + " psi=" + psi);
//			throw new Exception("Numerical error: p0 is non-positive " + top + "/" + bottom + "=" + result);
//		}
//		return result;
	}

	// Equation 2 from
	// Stadler, Tanja. "Sampling-through-time in birth–death trees." Journal of theoretical biology 267.3 (2010): 396-404.
	private double getP1(double height, double rho, double c1, double c2) throws Exception {
		return (4 * rho) / (2 * (1 - c2 * c2) + Math.exp(-c1 * height) * (1 - c2) * (1 - c2) + Math.exp(c1 * height) * (1 + c2) * (1 + c2));
//		double top = 4*rho;
//		double bottom = 2*(1-c2*c2) + Math.exp(-c1*height)*(1-c2)*(1-c2) + Math.exp(c1*height)*(1+c2)*(1+c2);
//		double result = top/bottom;
//		if (Double.isNaN(result) || result <= 0) {
//			//Log.warning("c1=" + c1 + " c2=" + c2 + " lambda =" + lambda + " mu=" + mu + " psi=" + psi);
//			throw new Exception("Numerical error: p1 is non-positive " + top + "/" + bottom + "=" + result);
//		}
//		return result;
	}


	// Equation 7 from
	// Douglas et al. "Evolution is coupled with branching across many granularities of life" (2025)
	// The (log) probability density of observing B stubs at ages Z along a single branch branchNr
	private double getStubLogPForBranchWithSampling(int branchNr, double lambda, double mu, double psi, double rho) {
		double stubLogP = 0;
		Tree tree = getTree();
		Stubs stubs = stubsInput.get();
		double treeHeight = tree.getRoot().getHeight(); // Height of the tree root
		Node node = tree.getNode(branchNr); // The node corresponding to branch branchNr

		// Check whether integrating over stubs (stub-free inference mode) is used
		if (!stubs.estimateStubs()) return 0;

		// Get the origin time h0 and end time h1 of branch branchNr
		double h0 = node.isRoot() ? treeHeight : node.getParent().getHeight();
		double h1 = node.getHeight();

		// If reversible jump is off
		// The number of stubs per branch is estimated
		if (!stubs.getReversibleJump()) {
			if (node.isRoot()) return 0; // If root node, return zero
			int stubCount = stubs.getNStubsOnBranch(branchNr);
			if (h0 <= h1) {
				if (stubCount == 0) return 0; // A zero-length branch must have zero stubs
				return Double.NEGATIVE_INFINITY; // Otherwise return negative infinity (impossible)
			}

			double expectedStubCount = getExpectedStubCountForBranch(h1, h0, lambda, mu, psi, rho); // Expected number of unobserved bifurcations over branchNr, E(B)

			// Poisson(expectedStubCount) distribution
			stubLogP = stubCount * Math.log(expectedStubCount) - expectedStubCount;
			for (int i = 2; i <= stubCount; i ++) stubLogP += -Math.log(i);

			return stubLogP;

		}

		// Else (if reversible jump is enabled)
		// The placement of each stub is estimated
		int stubCount = 0;
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

				double logRate = this.getStubLogBirthRate(stubHeight, lambda, mu, psi, rho);
				stubLogP += logRate;
				stubCount ++;
			}
		}

		// Integral and combinatorial term across the branch
		double expectedStubCount = getExpectedStubCountForBranch(h1, h0, lambda, mu, psi, rho);
		//Log.warning("g=" + g);
		stubLogP += -expectedStubCount;
		for (int i = 2; i <= stubCount; i ++) stubLogP += -Math.log(i);

		return stubLogP;
	}


	// Stub probability calculation when psi = 0 and rho = 1
	private double getStubLogPForBranchWithoutSampling(int branchNr, double lambda, double mu) {
		double stubLogP = 0;
		Tree tree = getTree();
		Stubs stubs = stubsInput.get();
		double treeHeight = tree.getRoot().getHeight(); // Height of the tree root
		Node node = tree.getNode(branchNr); // The node corresponding to branch branchNr

		// Check whether integrating over stubs (stub-free inference mode) is used
		if (!stubs.estimateStubs()) return 0;

		// Get the origin time h0 and end time h1 of branch branchNr
		double h0 = node.isRoot() ? treeHeight : node.getParent().getHeight();
		double h1 = node.getHeight();

		// If reversible jump is off
		// The number of stubs per branch is estimated
		if (!stubs.getReversibleJump()) {
			if (node.isRoot()) return 0; // If root node, return zero
			int stubCount = stubs.getNStubsOnBranch(branchNr);
			if (h0 <= h1) {
				if (stubCount == 0) return 0; // A zero-length branch must have zero stubs
				return Double.NEGATIVE_INFINITY; // Otherwise return negative infinity (impossible)
			}

			double expectedStubCount = getExpectedStubCountForBranch(h1, h0, lambda, mu);

			// Poisson(expectedStubCount) distribution
			stubLogP = stubCount * Math.log(expectedStubCount) - expectedStubCount;
			for (int i = 2; i <= stubCount; i ++) stubLogP += -Math.log(i);

			return stubLogP;

		}

		// Else (if reversible jump is enabled)
		// The placement of each stub is estimated
		int stubCount = 0;
		for (int stubNr = 0; stubNr < stubs.getStubDimension(); stubNr++) {
			if (!stubs.includeStub(stubNr)) continue;
			if (stubs.getBranches().getValue(stubNr) == branchNr) {
				double stubHeight = stubs.getAbsoluteTimeOfStub(stubNr);
				//double stubHeight = stubs.getRelativeTimeOfStub(stubNr);

				if (stubHeight < node.getHeight() || stubHeight > node.getParent().getHeight()) {
					return Double.NEGATIVE_INFINITY;
				}

				double logRate = this.getStubLogBirthRate(stubHeight, lambda, mu);
				stubLogP += logRate;
				stubCount ++;
			}
		}

		// Integral and combinatorial term across the branch
		double expectedStubCount = getExpectedStubCountForBranch(h1, h0, lambda, mu);
		stubLogP += -expectedStubCount;
		for (int i = 2; i <= stubCount; i ++) stubLogP += -Math.log(i);

		return stubLogP;

	}


	private double origin() {
		return originInput.get().getValue();
	}

	public double getLambda() {
		if (this.lambdaInput.get() == null) {
			double d = netDiversificationRateInput.get().getArrayValue();
			double r = turnoverInput.get().getArrayValue();
			return d / (1 - r);

			// Calculate without calling getMu() or there will be a never ending loop
//			if (r0Input.get() != null) {
//				double r = r0Input.get().getArrayValue();
//				return net / (1 - 1/r);
//			} else {
//				double t = turnoverInput.get().getArrayValue();
//				return net / (1 - t);
//			}

		} else {
			return this.lambdaInput.get().getArrayValue();
		}
	}

	public double getMu() {
		if (this.r0Input.get() == null) {
			double lambda = this.getLambda();
			double t = turnoverInput.get().getArrayValue();
			return t * lambda;
		} else {
			double lambda = this.getLambda();
			double r = r0Input.get().getArrayValue();
			return lambda / r;
		}
//		double lambda = this.getLambda();
//		if (r0Input.get() != null) {
//			return lambda / r0Input.get().getArrayValue();
//		} else {
//			return turnoverInput.get().getArrayValue() * lambda;
//		}
	}

	public double getPsi() {
		
		// Default
		if (psiInput.get() == null && samplingProportionInput.get() == null) {
			return 0.0;
		}
		
		if (psiInput.get() == null) {
			double mu = this.getMu();
			double s = samplingProportionInput.get().getArrayValue();
			return s * mu / (1 - s);
		} else {
			return this.psiInput.get().getArrayValue();
		}
//		if (samplingProportionInput.get() == null) return 0.0;
//		double mu = this.getMu();
//		double samplingProportion = samplingProportionInput.get().getArrayValue();
//		return samplingProportion*mu / (1-samplingProportion);
	}

	public double getRho() {
		
		// Default
		if (rhoInput.get() == null) {
			return 1.0;
		}
		return rhoInput.get().getValue();
	}

	public Tree getTree() {
		return (Tree) treeInput.get();
	}


	// Equation 1 from
	// Stolz, Stadler & Vaughan (2024) "Integrating Transmission Dynamics and Pathogen Evolution Through a Bayesian Approach" bioRxiv
	// https://doi.org/10.1101/2024.04.15.589468
	// The expected number of unobserved bifurcations on a branch between its origin time h0 and end time h1
	// when psi > 0 and rho < 1
	private double getExpectedStubCountForBranch(double h1, double h0, double lambda, double mu, double psi, double rho) {
		double c1 = getC1(lambda, mu, psi);
		double c2 = getC2(lambda, mu, psi, rho);
		double f = ((c2 - 1) * Math.exp(-c1 * h1) - c2 - 1) / ((c2 - 1) * Math.exp(-c1 * h0) - c2 - 1);
		return (h0 - h1) * (lambda + mu + psi - c1) + 2 * Math.log(f);
	}

	// The expected number of unobserved bifurcations on a branch between its origin time h0 and end time h1
	// when psi = 0 and rho = 1
	private double getExpectedStubCountForBranch(double h1, double h0, double lambda, double mu) {
		Tree tree = getTree();
		double treeHeight = tree.getRoot().getHeight();
		return getStubRateIntegral(h1, lambda, mu, treeHeight) - getStubRateIntegral(h0, lambda, mu, treeHeight);
	}
	// The following class used by getExpectedStubCountForBranch when psi = 0 and rho = 1
	private double getStubRateIntegral(double height, double lambda, double mu, double rootHeight) {
		double t = rootHeight - height; // Time elapsed since the root
		double b = Math.log(lambda * Math.exp((lambda - mu) * height) - mu);
//		return getStubMultiplier() * 2 * (b + lambda * s);
		return 2 * (b + lambda * t);
	}


//	/**
//	 * Multiply stub rate by this term
//	 * @return
//	 */
//	public double getStubMultiplier() {
//		// When no input for p, multiplier = 1
//		if (pInput.get() == null) return 1.0;
//		if (pInput.get().getDimension() == 1 && !perBranchSpikeInput.get()) {
//			return 1-pInput.get().getValue()/2;
//		} else if (pInput.get().getDimension() == 1 && perBranchSpikeInput.get()) {
//			return pInput.get().getValue();
//		} else {
//			double p1 = pInput.get().getValue(0); // Single mutant
//			double p2 = pInput.get().getValue(1); // Double mutant
//			double p3 = pInput.get().getValue(2); // No mutant
//			return p1/2 + p2;
//		}
//	}


	// The following two classes are used when the placement of stubs is estimated (reversible jump is enabled)
	// Equivalent to 2 * lambda * p0(t)
	// "Stubs appear in a time-dependent rate 2 * lambda * p0(t) in T."
	// psi > 0 and rho < 1
	private double getStubLogBirthRate(double timeRemaining, double lambda, double mu, double psi, double rho) {
		double logQt = this.getLogQt(timeRemaining, lambda, mu, psi, rho);
//		return Math.log(getStubMultiplier()) + Math.log(2) + Math.log(lambda) + logQt;
		return Math.log(2) + Math.log(lambda) + logQt; // 2 * lambda * p0(t)
	}
	// psi = 0 and rho = 1
	private double getStubLogBirthRate(double timeRemaining, double lambda, double mu) {
		double logQt = this.getLogQt(timeRemaining, lambda, mu);
//		return Math.log(getStubMultiplier()) + Math.log(2) + Math.log(lambda) + logQt;
		return Math.log(2) + Math.log(lambda) + logQt; // 2 * lambda * p0(t)
	}


	// Equation 1 from
	// Stadler, Tanja. "Sampling-through-time in birth–death trees." Journal of theoretical biology 267.3 (2010): 396-404.
	// Equivalent to p0(t)
	// when psi =! 0
	private double getLogQt(double time, double lambda, double mu, double psi, double rho) {
		double c1 = getC1(lambda, mu, psi);
		double c2 = getC2(lambda, mu, psi, rho);
		double exp = Math.exp(-c1 * time);
		double top = exp*(1-c2) - (1+c2);
		double btm = exp*(1-c2) + (1+c2);
		return Math.log(lambda + mu + psi + c1 * top/btm) - Math.log(2 * lambda);
	}

	// See Remark 3.2 from
	// Stadler, Tanja. "Sampling-through-time in birth–death trees." Journal of theoretical biology 267.3 (2010): 396-404.
	// Equivalent to p0(t)
	// when psi = 0 and rho = 1
	private double getLogQt(double time, double lambda, double mu) {
		double exp = Math.exp((lambda - mu) * time);
		double top = mu * (exp - 1);
		double btm = lambda * exp - mu;
        return Math.log(top) - Math.log(btm);
	}



	@Override
    protected boolean requiresRecalculation() {

//		if (pInput.get() != null && pInput.get().somethingIsDirty()) {
//			return true;
//		}

		return true;
//        return super.requiresRecalculation() ||
//        		InputUtil.isDirty(lambdaInput) ||
//				InputUtil.isDirty(muInput) ||
//				InputUtil.isDirty(psiInput) ||
//        		InputUtil.isDirty(r0Input) ||
//        		InputUtil.isDirty(netDiversificationRateInput) ||
//        		InputUtil.isDirty(turnoverInput) ||
//        		InputUtil.isDirty(samplingProportionInput) ||
//				InputUtil.isDirty(rhoInput) ||
//        		InputUtil.isDirty(stubsInput) ||
//				InputUtil.isDirty(originInput);
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



	@Override
	public double getMeanStubNumber(double nodeHeight, double parentHeight) {

		if (parentHeight <= nodeHeight) return 0;

		double lambda = this.getLambda();
		double mu = this.getMu();
		double psi = this.getPsi();
		double rho = this.getRho();

		// Calculate the expected number of unobserved bifurcations for a branch
		if (psi == 0 && rho == 1) {
			return getExpectedStubCountForBranch(nodeHeight, parentHeight, lambda, mu);
		} else {
			return getExpectedStubCountForBranch(nodeHeight, parentHeight, lambda, mu, psi, rho);
		}

	}


}