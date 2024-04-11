package gammaspike.tree;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

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
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;


@Description("A branch with length zero on a tree")
public class Stubs extends CalculationNode implements Loggable, Function {

	final public Input<TreeInterface> treeInput = new Input<>("tree", "tree without stubs", Input.Validate.REQUIRED);
	final public Input<IntegerParameter> branchNrInput = new Input<>("branchNr", "node index of each stub", Input.Validate.REQUIRED);
	final public Input<RealParameter> timeInput = new Input<>("time", "relative time of each stub along its branch", Input.Validate.REQUIRED);
	final public Input<IntegerParameter> labelIndicatorInput = new Input<>("indicator", "label indicator (-1, 0, 1)", Input.Validate.OPTIONAL);
	
	
	
	
	final public Input<RealParameter> originInput = new Input<>("origin", "length of origin branch", Validate.OPTIONAL);
	
	
	// Stochastic variable selection
	final public Input<BooleanParameter> selectionInput = new Input<>("selection", "boolean indicating if each stub is being used (stochastic variable selection)", Input.Validate.OPTIONAL);
	final public Input<Integer> maxNrStubsInput = new Input<>("maxNrStubs", "if set, then will use stochastic variable selection rather than reversible jump", -1);
	
	
	
	// For restoring
	double[] storedTimes;
	Integer[] storedBranches;
	Integer[] storedLabelIndicators;
	Boolean[] storedSelections;
	boolean reversibleJump;
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
		branchNrInput.get().setLower(0);
		branchNrInput.get().setUpper(nNodes-2); // Exclude the root
		branchNrInput.get().setValue(0);
		timeInput.get().setLower(0.0);
		timeInput.get().setUpper(1.0);
		timeInput.get().setValue(0.5);
		
		
		IntegerParameter labelIndicator = labelIndicatorInput.get();
		if (labelIndicator != null) {
			labelIndicator.setLower(-1);
			labelIndicator.setUpper(1);
			labelIndicator.setValue(0);
		}
		
		
		
		// Stochastic variable selection - fixed num stubs, but some are inactive
		if (maxNrStubsInput.get() > 0) {
			reversibleJump = false;
			
			//Log.warning(maxNrStubsInput.get() + " | " + selectionInput.get());
			if (selectionInput.get() == null) throw new IllegalArgumentException("Please provide 'selection' or set maxNrStubs to negative value");
			
			// Index 0 is the dummy
			int dim = maxNrStubsInput.get();
			branchNrInput.get().setDimension(dim+1);
			timeInput.get().setDimension(dim+1);
			selectionInput.get().setDimension(dim+1);
			
			// Initialise all at false 
			for (int i = 0; i < dim+1; i++) {
				selectionInput.get().setValue(i, false);
				//selectionInput.get().setValue(i, true);
			}
			
		}
		
		// Reversible jump - variable num stubs
		else {
			reversibleJump = true;
			
			// Initalise all vectors at length 1. This is the dummy position to ensure the lists are not empty
			branchNrInput.get().setDimension(1);
			timeInput.get().setDimension(1);
			if (labelIndicator != null) {
				labelIndicator.setDimension(1);
			}
		}

		
		
	}
	
	
	/**
	 * Return immutable list of stubs
	 * @return
	 */
	public List<Stub> getStubs(){
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
		if (this.reversibleJump) {
			return true;
		}
		return selectionInput.get().getValue(stubNr);
	}

	@Override
	public void init(PrintStream out) {
		out.print("nstubs\tstub1Height\t");
	}

	@Override
	public void log(long sample, PrintStream out) {
		out.print(this.getNStubs() + "\t");
		
		// Return the first stub time only
		double time = this.getNStubs() == 0 ? 0 : getAbsoluteTimeOfStub(1);
		//double time = this.getNStubs() == 0 ? 0 : getRelativeTimeOfStub(1);
		out.print(time + "\t");
		
		
	}

	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
		
	}
	
	public int getNStubs() {
		if (this.reversibleJump) return branchNrInput.get().getDimension() - 1; // Remove 1 for dummy index
		int nstubs = 0;
		for (int i = 1; i < branchNrInput.get().getDimension(); i ++) {
			if (this.includeStub(i)) nstubs ++;
		}
		return nstubs;
	}
	
	
	


	public int getNStubsOnBranch(int nodeNr) {
		int nstubs = 0;
		for (int i = 0; i < branchNrInput.get().getDimension(); i ++) {
			int b = branchNrInput.get().getValue(i);
			if (this.includeStub(i) && nodeNr == b) nstubs ++;
		}
		return nstubs;
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
		return this.reversibleJump;
	}
	
	
	public int getStubDimension() {
		return branchNrInput.get().getDimension();
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


	

}









