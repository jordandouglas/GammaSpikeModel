package gammaspike.tree;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import beast.base.core.BEASTInterface;
import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.CalculationNode;
import beast.base.inference.Distribution;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import gammaspike.distribution.BranchSpikePrior;
import gammaspike.distribution.StubExpectation;
import gammaspike.distribution.StumpedTreePrior;


@Description("A branch with length zero on a tree")
public class Stubs extends CalculationNode implements Loggable, Function {

	
	public enum StubMode{
		exact, 
		branchCounts,
		nostub
	}
	
	final public Input<TreeInterface> treeInput = new Input<>("tree", "tree without stubs", Input.Validate.REQUIRED);
	
	
	// Reversible jump MCMC
	final public Input<IntegerParameter> branchNrInput = new Input<>("branchNr", "node index of each stub", Input.Validate.OPTIONAL);
	final public Input<RealParameter> timeInput = new Input<>("time", "relative time of each stub along its branch", Input.Validate.OPTIONAL);
	final public Input<IntegerParameter> labelIndicatorInput = new Input<>("indicator", "label indicator (-1, 0, 1)", Input.Validate.OPTIONAL);
	
	
	// Integrate over stub heights
	final public Input<IntegerParameter> stubsPerBranchInput = new Input<>("stubsPerBranch", "number of stubs per branch", Input.Validate.OPTIONAL);
	
	
	final public Input<StubExpectation> priorInput = new Input<>("prior", "the prior object for calculating stub expectations", Validate.OPTIONAL);
	
	final public Input<RealParameter> originInput = new Input<>("origin", "length of origin branch", Validate.OPTIONAL);
	

	StubExpectation stubExpectation = null;
	BranchSpikePrior spikePrior = null;
	StubMode stubMode;
	long[] sampledStubSampleNr = null;
	int[] sampledStubNr = null;
	
	// For restoring
	double[] storedTimes;
	Integer[] storedBranches;
	Integer[] storedLabelIndicators;
	Boolean[] storedSelections;
	//Integer[] storedTos;
	//Integer[] storedFroms;
		
	
	public class Stub {
		
		private int stubIndex;
		private int branchNr;
		private double absoluteHeight;
		private int labelIndicator;
		
		public Stub(int stubIndex, int branchNr, double absoluteHeight, int labelIndicator) {
			this.stubIndex = stubIndex;
			this.branchNr = branchNr;
			this.absoluteHeight = absoluteHeight;
			this.labelIndicator = labelIndicator;
		}
		
		public int getBranchNr() {
			return this.branchNr;
		}
		
		public double getAbsoluteHeight() {
			return this.absoluteHeight;
		}
		
		
		/**
		 * -1 = stub is ancestral
		 *  0 = subfunctionalisation
		 * +1 = stub is mutant
		 * @return
		 */
		public int getLabelIndicator() {
			return this.labelIndicator;
		}
		
		
		@Override
		public String toString() {
			return "Stub at " + this.getAbsoluteHeight() + " on node " + branchNr;
		}
		
		
		// String for tree metadata
		public void toMetaData(StringBuffer buf, int index) {
			buf.append(getStubTimeName(index)).append("=").append(absoluteHeight).append(",");
			buf.append(getStubFromSymbolName(index)).append("='").append(labelIndicator).append("',");
			buf.append(getStubToSymbolName(index)).append("='").append("").append("'");
		}

		public int getIndex() {
			return stubIndex;
		}

	}
	

	@Override
	public void initAndValidate() {
		
		
		
		
		
		// Lower and uppers
		int nNodes = treeInput.get().getNodeCount();

		
		IntegerParameter labelIndicator = labelIndicatorInput.get();
		if (labelIndicator != null) {
			labelIndicator.setLower(-1);
			labelIndicator.setUpper(1);
			labelIndicator.setValue(0);
		}
		
		
		if (priorInput.get() == null && stubsPerBranchInput.get() == null && branchNrInput.get() == null) {
			throw new IllegalArgumentException("Please provide either prior, stubsPerBranch, or branchNr+time");
		}
		
		// No stubs - they are integrated over
		if (priorInput.get() != null && stubsPerBranchInput.get() == null && branchNrInput.get() == null) {
		
			this.stubMode = StubMode.nostub;
			
			this.sampledStubNr = new int[nNodes];
			this.sampledStubSampleNr = new long[nNodes];
			for (int i = 0; i < sampledStubSampleNr.length; i ++) {
				sampledStubSampleNr[i] = -1;
			}
			
			
			stubExpectation = priorInput.get();
			Log.warning("Found distribution " + ((BEASTObject)stubExpectation).getID());
			
			// Find the object that is used to compute the spike prior
			for (BEASTInterface o : this.getOutputs()) {
				if (o != null && o instanceof BEASTObject) {
					//Log.warning("input :" + o.getClass().getCanonicalName());
					
					if (o instanceof BranchSpikePrior) {
						
						if (spikePrior != null) {
							throw new IllegalArgumentException("Found multiple BranchSpikePrior objects that are pointing to "
									+ this.getID() + ". Please ensure there is 1 at most (or estimate stubs directly).");
						}
						
						spikePrior = (BranchSpikePrior)o;
						Log.warning("Found spike prior " + spikePrior.getID());
					}
					
					else if (o instanceof Distribution && o != stubExpectation) {
						throw new IllegalArgumentException("Found an unknown distribtuion " + o.getID() + " pointing to "
								+ this.getID() + ". Please ensure the only priors being used are the tree prior and a spike prior (or estimate stubs directly).");
				
						
					}
				}
				
				
			}
			
			if (stubExpectation == null) {
				throw new IllegalArgumentException("Cannot find an object that implements 'StubExpectation'");
			}
			
		
		}
		
		// The stub count on each branch is estimated
		else if (stubsPerBranchInput.get() != null) {
			
			this.stubMode = StubMode.branchCounts;
			
			IntegerParameter p = stubsPerBranchInput.get();
			
			// Do not initialise if we are resuming
			p.setLower(0);
			p.setDimension(nNodes-2);
			for (int i = 0; i < p.getDimension(); i ++) {
				p.setValue(i, 0);
			}
			
		}
		
		// The branch and height of each stis estimated
		else {
			this.stubMode = StubMode.exact;
			
			boolean reset = true;
			
			branchNrInput.get().setLower(0);
			branchNrInput.get().setUpper(nNodes-2); // Exclude the root
			branchNrInput.get().setValue(0);
			timeInput.get().setLower(0.0);
			timeInput.get().setUpper(1.0);
			timeInput.get().setValue(0.5);
			
			
			// Do not initialise if we are resuming
			if (timeInput.get().getDimension() > 1 & 
					branchNrInput.get().getDimension() > 1 & 
					timeInput.get().getDimension() == timeInput.get().getDimension()) {
				
				reset = false;
				
				Log.warning("Resuming stubs " + getNStubs());
				
			}
			
			if (reset) {
				
				
				Log.warning("Initialising stubs" + getNStubs() + timeInput.get().getDimension() + " . " + branchNrInput.get().getDimension());
			
				// Initalise all vectors at length 1. This is the dummy position to ensure the lists are not empty
				branchNrInput.get().setDimension(1);
				timeInput.get().setDimension(1);
				if (labelIndicator != null) {
					labelIndicator.setDimension(1);
				}
			
			}
		}

		
		
	}
	
	
	/**
	 * Return immutable list of stubs
	 * @return
	 */
	public List<Stub> getStubs(){
		
		if (this.stubMode == StubMode.nostub) {
			throw new IllegalArgumentException("Dev error: no stubs but getStubs() was called");
		}
		
		List<Stub> stubs = new ArrayList<>();
		for (int i = 1; i <= this.getStubDimension(); i ++) {
			int branchNr = this.getBranches().getValue(i);
			double time = getAbsoluteTimeOfStub(i);
			boolean include = this.includeStub(i);
			int labelIndicator = !hasIndicatorLabels() ? 0 : labelIndicatorInput.get().getValue(i);
			if (include) {
				Stub stub = new Stub(i, branchNr, time, labelIndicator);
				stubs.add(stub);
			}
		}
		return stubs;
	}
	
	
	/**
	 * Stub is always included in reversible jump, and sometimes included in SVS
	 * @param stubNr
	 * @return
	 */
	public boolean includeStub(int stubNr) {
		if (stubNr == 0) return false;
		return true;
		//if (this.reversibleJump) {
			//return true;
		//}
		//return selectionInput.get().getValue(stubNr);
	}

	@Override
	public void init(PrintStream out) {
		out.print("nstubs\t");
	}

	@Override
	public void log(long sample, PrintStream out) {
		if (!this.estimateStubs()) {
			out.print(this.sampleNStubs(sample) + "\t");
		}else {
			out.print(this.getNStubs() + "\t");
		}
		
	}

	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
		
	}
	
	public int getNStubs() {
		
		
		switch (this.stubMode) {
			case nostub:
				
				return -1;
				
			case branchCounts:
				
				int nstubs = 0;
				for (int i = 0; i < this.stubsPerBranchInput.get().getDimension(); i++) {
					nstubs += stubsPerBranchInput.get().getValue(i);
				}
				return nstubs;
				
			case exact:
				
				return branchNrInput.get().getDimension() - 1; // Remove 1 for dummy index
				
		}
		
		return -1;

	}
	
	
	/**
	 * Sample all nodes on the tree, for logging
	 * @param sampleNr
	 * @return
	 */
	public int sampleNStubs(long sampleNr) {
		
		if (this.estimateStubs()) return -1;
		
		int nstubs = 0;
		for (int i = 0; i < this.treeInput.get().getNodeCount(); i++) {
			nstubs += sampleNStubsOnBranch(i, sampleNr);
		}
		
		
		//Log.warning(" sampleNStubs " + sampleNr + " " + nstubs);
		return nstubs;
		
		
	}
	
	
	/**
	 * Sample the number of nodes on this branch, for logging
	 * @param nodeNr
	 * @param sampleNr
	 * @return
	 */
	public int sampleNStubsOnBranch(int nodeNr, long sampleNr) {
		
		if (this.estimateStubs()) return -1;
		
		
		// No stubs on root
		if (nodeNr == treeInput.get().getRoot().getNr()) {
			return 0;
		}
		
		// To enssyre harmony across all loggers, don't sample again on this state if it has already been sampled
		if (sampleNr == this.sampledStubSampleNr[nodeNr]) {
			return this.sampledStubNr[nodeNr];
		}
		
		
		
		
		Node node = treeInput.get().getNode(nodeNr);
		double h0 = node.getHeight();
		double h1 = node.isRoot() ? node.getHeight() : node.getParent().getHeight(); 
		double poissonMean = stubExpectation.getMeanStubNumber(h0, h1);
		
		int nstubs = 0;
		
		if (poissonMean == 0) {
			nstubs = 0;
		}
		
		else if (this.spikePrior == null) {
			
			// Sample from Poisson distribution if there are no spikes
			nstubs = (int) Randomizer.nextPoisson(poissonMean);
			
			
		}else {
			
			// Sample from a Gamma-Poisson distribution
			double[] cf = this.spikePrior.getCumulativeProbs(poissonMean, nodeNr);
			nstubs = Randomizer.randomChoice(cf);
			
		}
		
		
		//Log.warning(" on branch " + h1 + ", " + node.getHeight() + " mean is " + poissonMean + " and n is " + nstubs);
		
		
		this.sampledStubSampleNr[nodeNr] = sampleNr;
		this.sampledStubNr[nodeNr] = nstubs;
		return nstubs;
		
		
	}
	

	public int getNStubsOnBranch(int nodeNr) {
		
		
		switch (this.stubMode) {
		
		// Stochastic function
		case nostub:
			return -1;
			
		case branchCounts:
			
			if (nodeNr >= this.stubsPerBranchInput.get().getDimension()) return 0;
			return this.stubsPerBranchInput.get().getValue(nodeNr);
			
		case exact:
			
			int nstubs = 0;
			for (int i = 0; i < branchNrInput.get().getDimension(); i ++) {
				int b = branchNrInput.get().getValue(i);
				if (this.includeStub(i) && nodeNr == b) nstubs ++;
			}
			return nstubs;
			
		}
	
		return -1;

		
	}
	

	public int getNStubsOnBranchWithIndicator(int nodeNr, int indicatorValue) {
		int nstubs = 0;
		for (int i = 0; i < branchNrInput.get().getDimension(); i ++) {
			int b = getBranches().getValue(i);
			int ind = getLabelIndicators().getValue(i);
			if (this.includeStub(i) && nodeNr == b && ind == indicatorValue) nstubs ++;
		}
		return nstubs;
	}
	
	
	
	public IntegerParameter getBranches() {
		return branchNrInput.get();
	}

	
	public RealParameter getStubHeights() {
		return timeInput.get();
	}
	
	public IntegerParameter getLabelIndicators() {
		return labelIndicatorInput.get();
	}
	
	/**
	 * Call by operator when dimensions change, since the default store/restore does not account for dimensional change
	 */
	public void storeDimensions() {
		storedTimes = getStubHeights().getDoubleValues();
		storedBranches = getBranches().getValues();
		if (getLabelIndicators() != null) {
			storedLabelIndicators = getLabelIndicators().getValues();
		}
		//if (this.selectionInput.get() != null) {
			//storedSelections = getSelection().getValues();
		//}
		
	}
	


	public void acceptDimensions() {
		storedTimes = null;
		storedBranches = null;
		storedLabelIndicators = null;
		//storedSelections = null;
		//super.accept();
	}
	
	

    /**
     * CalculationNode methods *
     */
    @Override
    public void store() {
        super.store();
    }

    @Override
    public void restore() {
    	restoreDimensions();
        super.restore();
        
    }
	
	private void restoreDimensions() {
		
		
		if (storedTimes == null) return;
		
		RealParameter times = this.getStubHeights();
		IntegerParameter branches = this.getBranches();
		IntegerParameter labelIndicators = this.getLabelIndicators();
		
		int oldDimension = storedTimes.length;
		
		
		//Log.warning("restore from " + times.getDimension() + " to "  + oldDimension);
		
		times.setDimension(oldDimension);
		branches.setDimension(oldDimension);
		if (labelIndicators != null) labelIndicators.setDimension(oldDimension);
		
		for (int i = 0; i < oldDimension; i ++) {
			times.setValue(i, storedTimes[i]);
			branches.setValue(i, storedBranches[i]);
			if (labelIndicators != null) labelIndicators.setValue(i, storedLabelIndicators[i]);
		}
		
		
		storedTimes = null;
		storedBranches = null;
		storedLabelIndicators = null;
		
	}



	/*
	 * Returns a list of concert event indices along this branch, sorted by time (forward in time)
	 */
	public List<Stub> getSortedStubsOnBranch(Node node){
		
		
		if (this.stubMode != StubMode.exact) {
			throw new IllegalArgumentException("DevError 42574: cannot sort stubs unless RJ is being used");
		}
		
		// Get indices and times. Start at index 1 because 0 is the dummy/placeholder
		List<Stub> stubs = new ArrayList<>();
		for (int stubNr = 1; stubNr < this.getStubDimension(); stubNr++) {
			
			if (!this.includeStub(stubNr)) continue;
			int branchNr = branchNrInput.get().getNativeValue(stubNr);
			
			if (branchNr == node.getNr()) {
				double time = getAbsoluteTimeOfStub(stubNr);
				int labelIndicator = !hasIndicatorLabels() ? 0 : labelIndicatorInput.get().getValue(stubNr);
				stubs.add(new Stub(stubNr, branchNr, time, labelIndicator));
			}
		}
		
		
		this.sortStubs(stubs);
		
		return stubs;
	}
	
	/**
	 * In place sorting of events by time (forward in time)
	 * @param events
	 */
	public void sortStubs(List<Stub> events) {
		
		// Sort the list by time (forward in time; ie. reverse order of height)
		Collections.sort(events, new Comparator<Stub>() {
		    public int compare(Stub left, Stub right) {
		        return Double.compare(right.getAbsoluteHeight(), left.getAbsoluteHeight());
		    }
		});
		
	}
	
	
	
	
	
	
	public static String getStubCountName() {
		return "ConcertCount";
	}

	public static String getStubTimeName(int eventNr) {
		return "ceTime" + eventNr;
	}
	
	public static String getStubFromSymbolName(int eventNr) {
		return "ceFrom" + eventNr;
	}
	
	public static String getStubToSymbolName(int eventNr) {
		return "ceTo" + eventNr;
	}
	
	public double getLabelIndicatorOfStub(int stubNr) {
		if (this.labelIndicatorInput.get() == null) return 0;
		return this.labelIndicatorInput.get().getValue(stubNr);
	}

	
	public boolean hasIndicatorLabels() {
		return this.labelIndicatorInput.get() != null;
	}
	
	

	public double getAbsoluteTimeOfStub(int stubNr) {
		
		// Convert relative time to true time
		int branchNr = getBranches().getValue(stubNr);
		Node node = treeInput.get().getNode(branchNr);
		double relativeTime = timeInput.get().getArrayValue(stubNr);
		double time;
		if (node.isRoot()) {
			time = node.getHeight() + originInput.get().getValue();
			throw new IllegalArgumentException("Stubs on root node are not supported. Node nr " + branchNr);
		}else {
			time = node.getHeight() + relativeTime*(node.getParent().getHeight() - node.getHeight());
		}
		return time;
		
	}


	public boolean getReversibleJump() {
		return this.stubMode == StubMode.exact;
	}
	
	public boolean estimateStubs() {
		return this.stubMode != StubMode.nostub;
	}
	
	public int getStubDimension() {
		if (!getReversibleJump()) {
			return stubsPerBranchInput.get().getDimension();
		}else {
			return branchNrInput.get().getDimension();
		}
		
	}
	

	
	// Function overrides
	
	@Override
	public double getArrayValue(int dim) {
		return this.getNStubs();
	}


	@Override
	public int getDimension() {
		return 1;
	}


	public double getRelativeTimeOfStub(int stubNr) {
		double relativeTime = timeInput.get().getArrayValue(stubNr);
		return relativeTime;
	}


	public boolean stubIsAncestral(int stubNr) {
		if (!this.hasIndicatorLabels()) return true;
		int indicator = this.getLabelIndicators().getValue(stubNr);
		if (indicator == -1 || indicator == 0) {
			return true;
		}
		return false;
	}


	/**
	 * Return a list of branch numbers that have at least 1 stub
	 * @return
	 */
	public List<Integer> getBranchesWithStubs() {
		List<Integer> branches = new ArrayList<>();
		for (int stubNr = 0; stubNr < this.getStubDimension(); stubNr++) {
			if (!this.includeStub(stubNr)) continue;
			int branchNr = this.getBranches().getValue(stubNr);
			if (!branches.contains(branchNr)) branches.add(branchNr);
		}
		return branches;
	}


	/**
	 * Sample an active stub uniformly at random
	 * @return
	 */
	public Stub getRandomStub() {
		int nstubs = this.getNStubs();
		int stubPosSampled = Randomizer.nextInt(nstubs);
		int i = 0;
		for (int stubNr = 0; i < this.getStubDimension(); stubNr++) {
			
			if (!this.includeStub(stubNr)) continue;
			if (i == stubPosSampled) {
				int branchNr = this.getBranches().getValue(i);
				double time = getAbsoluteTimeOfStub(i);
				int labelIndicator = !hasIndicatorLabels() ? 0 : labelIndicatorInput.get().getValue(i);
				return new Stub(stubNr, branchNr, time, labelIndicator);
			}
			i++;
		}
		
		
		return null;
	}


	public boolean hasIndicators() {
		return this.labelIndicatorInput.get() != null;
	}



	
	
	
	/**
	 * Cache all branch lengths before a tree proposal
	 * Call this before getLogJacobian()
	 */
	public double[] prepareJacobian() {
		
		if (!this.getReversibleJump()) return null;
		
		Tree tree = (Tree)treeInput.get();
		double[] cachedBranchLengths = new double[tree.getNodeCount()];
		for (int branchNr = 0; branchNr < tree.getNodeCount(); branchNr++) {
			cachedBranchLengths[branchNr] = tree.getNode(branchNr).getLength();
		}
		return cachedBranchLengths;
	}
	
	
	/**
	 * Calculate Jacobian - make sure to call prepareJacobian before making a tree proposal
	 * Applicable only if tree topolgy or lengths change, but stub parameters are not changed
	 * This is because stub node heights are relative to branch length. So changing the length of 
	 * a branch that a stub is on requires knowledge of the ratio between the before and after to calculate Jacobian
	 */
	public double getLogJacobian(double[] cachedBranchLengths) {
		
		if (!this.getReversibleJump()) return 0;
		if (cachedBranchLengths == null) {
			throw new IllegalArgumentException("Developer error 111231: please call 'prepareJacobian' before 'getJacobian'");
		}
		
		
		double logJacobian = 0;
		Tree tree = (Tree)treeInput.get();
		for (int stubNr = 0; stubNr < this.getStubDimension(); stubNr++) {
			
			
			if (!this.includeStub(stubNr)) continue;
        	int stubBranchNr = this.getBranches().getValue(stubNr);
        	
        	double oldLength = cachedBranchLengths[stubBranchNr];
        	double newLength = tree.getNode(stubBranchNr).getLength();
        	if (oldLength != newLength) {
        		logJacobian += Math.log(newLength/oldLength);
        	}
			
		}
		
		
		return logJacobian;
	}


	public double getMeanNumberOfStubs(double h0, double h1) {
		if (this.stubExpectation == null) return -1;
		return this.stubExpectation.getMeanStubNumber(h0, h1);
	}




	

}









