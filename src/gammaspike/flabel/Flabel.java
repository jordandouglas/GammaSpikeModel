package gammaspike.flabel;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import gammaspike.tree.ForwardTimeSimulatorResub;
import gammaspike.tree.Stubs;


@Description("Tracks the internal node labels")
public class Flabel extends CalculationNode implements Function, Loggable {

	
	public enum BurstModel{
		BINARY,
		THREE_STATE,
		FOUR_STATE
	}
	

	public enum FlabelIndicator{
		LEFT_BURST(0),
		RIGHT_BURST(1),
		DOUBLE_BURST(2),
		NO_BURST(3);

        private final int value;

        FlabelIndicator(final int value) {
            this.value = value;
        }

        public int val() { 
        	return value; 
        }
        
        public static FlabelIndicator get(int value) {
        	if (value == 0) return LEFT_BURST;
        	if (value == 1) return RIGHT_BURST;
        	if (value == 2) return DOUBLE_BURST;
        	if (value == 3) return NO_BURST;
        	return null;
        	
        }
    }
	
	final public Input<TreeInterface> treeInput = new Input<>("tree", "tree without stubs", Input.Validate.REQUIRED);
	final public Input<IntegerParameter> leafLabelsInput = new Input<>("data", "leaf labels, one number per taxon. be careful with ordering", Input.Validate.OPTIONAL);
	final public Input<List<FlabelData>> dataInput = new Input<>("d", "more precise way of specifying data", new ArrayList<>());
	
	
	final public Input<IntegerParameter> indicatorsInput = new Input<>("indicators", "internal node indicators (0, 1, 2, 3)", Input.Validate.REQUIRED);
	final public Input<Stubs> stubsInput = new Input<>("stubs", "stubs for the tree", Input.Validate.OPTIONAL);
	
	final public Input<Boolean> simulatingInput = new Input<>("simulating", "set to true if simulating the stub bursts", false);
	final public Input<IntegerParameter> stubsPerBranchInput = new Input<>("nstubsPerBranch", "number of stubs per branch (for simulation)", Input.Validate.OPTIONAL);
	final public Input<RealParameter> pInput = new Input<>("p", "the probability of a single-mutant", Input.Validate.REQUIRED);
	final public Input<Boolean> perBranchSpikeInput = new Input<>("perBranchSpike", "if true, then the indicator will disregard labels", false);

	
	
	BurstModel burstModelType;
	int nextLabelNr = 0;
	int nleaves;
	int nnodes;
	int[] nodeLabels;
	Stubs stubs;
	
	
	int[] nburstsPerBranch;
	int nBurstStubsTotal;
	
	
	
	@Override
	public void initAndValidate() {
		
		
		if (simulatingInput.get() && leafLabelsInput.get() == null) {
			throw new IllegalArgumentException("Please specify data when simulating");
		}
		
		if ((leafLabelsInput.get() == null && dataInput.get().isEmpty())
				&& (leafLabelsInput.get() == null && !dataInput.get().isEmpty())
				&& (leafLabelsInput.get() != null && dataInput.get().isEmpty())) {
			throw new IllegalArgumentException("Please specify either data or flabelData");
		}
		

		
		// Initialise indicator dimension and range
		int pDim = pInput.get().getDimension();
		if (pDim == 1 && perBranchSpikeInput.get()) {
			this.burstModelType = BurstModel.BINARY;
			indicatorsInput.get().setLower(2);
			indicatorsInput.get().setUpper(3);
		}else if (pDim == 1 && !perBranchSpikeInput.get()) {
			this.burstModelType = BurstModel.THREE_STATE;
			indicatorsInput.get().setLower(0);
			indicatorsInput.get().setUpper(2);
		}else if (pDim == 3) {
			this.burstModelType = BurstModel.FOUR_STATE;
			indicatorsInput.get().setLower(0);
			indicatorsInput.get().setUpper(3);
		}else {
			throw new IllegalArgumentException("Please ensure p has dimension 1 or 3");
		}
		
		
		this.stubs = stubsInput.get();
		this.nnodes = treeInput.get().getNodeCount();
		this.nleaves = treeInput.get().getLeafNodeCount();
		
		if (!simulatingInput.get()) {
			
			if (leafLabelsInput.get() != null) {
				leafLabelsInput.get().setLower(0);
				leafLabelsInput.get().setUpper(nleaves-1);
			}
			
			if (leafLabelsInput.get() != null && leafLabelsInput.get().getDimension() != nleaves) {
				throw new IllegalArgumentException("Expected " + nleaves + " leaf datapoints but got " + leafLabelsInput.get().getDimension());
			}
			if (leafLabelsInput.get() == null && dataInput.get().size() != nleaves) {
				throw new IllegalArgumentException("Expected " + nleaves + " data datapoints but got " + dataInput.get().size());
			}
		}
		

		
		this.nodeLabels = new int[nnodes];
		initialiseIndicators(false);
		

		
		// Sample the number of bursts along each branch
		if (simulatingInput.get()) {
			nBurstStubsTotal = 0;
			nburstsPerBranch = new int[nnodes];
			for (int i = 0; i < nnodes; i ++) nburstsPerBranch[i] = 0;
			
			
			if (stubsPerBranchInput.get() == null || pInput.get() == null) {
				throw new IllegalArgumentException("Please provide nstubsPerBranch and p if simulating=true");
			}
			
			double pStubInducesBurst = 0;
			
			
			if (this.burstModelType == BurstModel.BINARY) {
				pStubInducesBurst = pInput.get().getValue();
				Log.warning("2 add stub wp " + pStubInducesBurst);
			}else if (this.burstModelType == BurstModel.THREE_STATE) {
				pStubInducesBurst = 1-pInput.get().getValue()/2;
				Log.warning("3 add stub wp " + pStubInducesBurst);
			}else if (this.burstModelType == BurstModel.FOUR_STATE) {
				pStubInducesBurst = pInput.get().getValue(0)/2 + pInput.get().getValue(1);
				Log.warning("4 add stub wp " + pStubInducesBurst);
			}
			
			
			stubsPerBranchInput.get().setDimension(nnodes);
			for (int i = 0; i < nnodes; i ++) {
				
				
				Node node = treeInput.get().getNode(i);
				
				
				// Bursts along this branch
				int nbursts = 0;
				
				// Bursts per stub
				int nstubsOnBranch = (int) ((double) node.getMetaData(ForwardTimeSimulatorResub.NSTUBS_STR));
				
				stubsPerBranchInput.get().setValue(i, nstubsOnBranch);
				for (int stubNr = 0; stubNr < nstubsOnBranch; stubNr++) {
					if (Randomizer.nextFloat() < pStubInducesBurst) {
						nbursts ++;
						nBurstStubsTotal++;
					}
				}
				
				//Log.warning("node " + i + " has " + nbursts + " " + node.toNewick());
				nburstsPerBranch[i] = nburstsPerBranch[i] + nbursts;
				
				// Give one to either child?
				if (!node.isLeaf()) {
					
					FlabelIndicator indicator = getIndicatorOfNode(node);
					
					if (indicator == FlabelIndicator.NO_BURST) continue;
					
					
					// Mutate right child
					if (indicator == FlabelIndicator.RIGHT_BURST || indicator == FlabelIndicator.DOUBLE_BURST) {
						int rightNr = node.getRight().getNr();
						nburstsPerBranch[rightNr] = nburstsPerBranch[rightNr] + 1;
						
					}
					
					// Mutate left child
					if (indicator == FlabelIndicator.LEFT_BURST || indicator == FlabelIndicator.DOUBLE_BURST) {
						int leftNr = node.getLeft().getNr();
						nburstsPerBranch[leftNr] = nburstsPerBranch[leftNr] + 1;
					}
					
					
				}
				
			}
			
			
			
			// Sample leaf labels
			update();
			leafLabelsInput.get().setDimension(nleaves);
			for (int i = 0; i < nleaves; i ++) {
				leafLabelsInput.get().setValue(i, this.nodeLabels[i]);
			}
			
			
			
		}else {
			nBurstStubsTotal = this.stubs.getNStubs();
		}
		
		resetLeafLabels();
		
				
		
	}
	
	public void resetLeafLabels() {
		
		// Set leaf labels
		for (int leafNr = 0; leafNr < nleaves; leafNr++) {
			this.nodeLabels[leafNr] = getLeafLabel(treeInput.get().getNode(leafNr));
		}
	}
	
	
	public int getLeafLabel(Node node) {
		
		if (!node.isLeaf()) return -1;
		
		if (leafLabelsInput.get() != null) {
			return leafLabelsInput.get().getValue(node.getNr());
		}
		
		for (FlabelData datapoint : dataInput.get()) {
			if (datapoint.getTaxon().equals(node.getID())) {
				return datapoint.getLabel();
			}
		}
		
		Log.warning("Cannot find " + node.getID());
		return -1;
		
	}
	
		public BurstModel getBurstModel() {
		return this.burstModelType;
	}
	
	
	
	/**
	 * Initialise internal node indicators. Care is needed when assigning burst vs no burst
	 */
	public void initialiseIndicators(boolean allowError) {
		
		
		int ninternal = nleaves-1;
		
		// Initialise indicators, as all -1, 0, or 1
		if (this.burstModelType == BurstModel.BINARY) {
			Log.warning("Initialising indicators to length " + (ninternal + nleaves));
			indicatorsInput.get().setDimension(ninternal + nleaves);
		}else {
			Log.warning("Initialising indicators to length " + ninternal);
			indicatorsInput.get().setDimension(ninternal);
		}
		
		
		// In this model, there is one indicator per branch, but other models have one per internal node
		if (this.burstModelType == BurstModel.BINARY) {
			
			double p = pInput.get() == null ? 0.5 : pInput.get().getValue();
			for (int i = 0; i < indicatorsInput.get().getDimension(); i++) {
				double u = Randomizer.nextDouble();
				if (u < p) indicatorsInput.get().setValue(i, FlabelIndicator.DOUBLE_BURST.val());
				else indicatorsInput.get().setValue(i, FlabelIndicator.NO_BURST.val());
			}
			
			
		}
		
		
		else if (this.burstModelType == BurstModel.THREE_STATE)  {
			double p = pInput.get() == null ? 0.6667 : pInput.get().getValue();
			for (int i = 0; i < indicatorsInput.get().getDimension(); i++) {
				double u = Randomizer.nextDouble();
				if (u < p/2) indicatorsInput.get().setValue(i, FlabelIndicator.LEFT_BURST.val());
				else if (u < p) indicatorsInput.get().setValue(i, FlabelIndicator.RIGHT_BURST.val());
				else indicatorsInput.get().setValue(i, FlabelIndicator.DOUBLE_BURST.val());
			}
		}
		
		
		else if (this.burstModelType == BurstModel.FOUR_STATE) {
			
			
			// Check probabilities are valid
			double p1 = pInput.get().getValue(0); // Single mutant
			double p2 = pInput.get().getValue(1); // Double mutant
			double p3 = pInput.get().getValue(2); // No mutant
			double sum = p1 + p2 + p3;
			
			if (Math.abs(1-sum) > 1e-6) {
				throw new IllegalArgumentException("Please ensure that p sums to 1, currently sums to " + sum);
			}
			
			
			if (simulatingInput.get()) {
				
				// Sample from dist
				for (int i = 0; i < indicatorsInput.get().getDimension(); i++) {
					double u = Randomizer.nextDouble()*sum;
					if (u < p1/2) indicatorsInput.get().setValue(i, FlabelIndicator.LEFT_BURST.val());
					else if (u < p1) indicatorsInput.get().setValue(i, FlabelIndicator.RIGHT_BURST.val());
					else if (u < p1+p2) indicatorsInput.get().setValue(i, FlabelIndicator.DOUBLE_BURST.val());
					else indicatorsInput.get().setValue(i, FlabelIndicator.NO_BURST.val());
					
					
					//Log.warning(u + " " + p1 + " " + p2 + " " + p3 + " " + indicatorsInput.get().getValue(i));
					
				}
					
				
			}else {
				
				
				final int nattempts = 1000;
				for (int attempt = 1; attempt <= nattempts; attempt++) {
					if (attempt >= nattempts) {
						if (allowError) {
							throw new IllegalArgumentException("Cannot initialise labelling after " + nattempts + " attempts");
						}else {
							return;
						}
					}
					try {
						resetMaxLabelNr();
						initialiseIndicators(treeInput.get().getRoot());
						if (isValid()) {
							break;
						}
					}catch(Exception e) {
						e.printStackTrace();
					}
				}
				
			}
			
			
		}
			
			
			
		
		
			
		
	}
	
	
	private void resetMaxLabelNr() {
		nextLabelNr = 0;
		for (int i = 0; i < nleaves; i++) {
			nextLabelNr = Math.max(nextLabelNr, getLeafLabel(treeInput.get().getNode(i)));
		}
		nextLabelNr++;
	}


	private int initialiseIndicators(Node node) throws Exception {
		
		if (node.isLeaf()) {
			int leafLabel = getLeafLabel(node);
			//Log.warning(node.getID() + " has " + leafLabel);
			return leafLabel;
		}

		
		// What are labels of children?
		int leftLabel = initialiseIndicators(node.getLeft());
		int rightLabel = initialiseIndicators(node.getRight());
		
		//Log.warning("LR=" + leftLabel + "," + rightLabel);
		
		
		// Are there any stubs on child branches?
		Stubs stubs = stubsInput.get();
		int nstubs = stubs.getNStubsOnBranch(node.getLeft().getNr()) + stubs.getNStubsOnBranch(node.getRight().getNr());

		
		int i = node.getNr() - this.nleaves;
		
		
		// Confirm that there does not exist a label in the leaf sets of both children, as we are assuming convergence is impossible
		boolean mutuallyExclusiveLabels = true;
		int commonLabel = -1;
		List<Integer> leftLabels = this.getLabelSet(node.getLeft());
		List<Integer> rightLabels = this.getLabelSet(node.getRight());
		for (int x : leftLabels) {
			if (rightLabels.contains(x)) {
				mutuallyExclusiveLabels = false;
				commonLabel = x;
				break;
			}
		}
		
		
		// No burst is necessary if the two labels are the same and there are no stubs
		if (!mutuallyExclusiveLabels || (leftLabel == rightLabel && nstubs == 0)) {
			indicatorsInput.get().setValue(i, FlabelIndicator.NO_BURST.val());
			
			if (leftLabel != commonLabel && rightLabel != commonLabel) {
				
				// Invalid
				//throw new Exception("Bad init");
				
				
			}
			
			//Log.warning("setting to nb " + leftLabel + " " + node.toNewick());
			return leftLabel; // = rightLabel
		}
		
		
		// There MAY have been a burst, or there MUST have been a burst, depending on the situation. But we will assign one either way
		else {
			
			
			// If the label on the left (right) is found elsewhere in the tree, then this node must be a burst right (left)
			Node parent = node.getParent();
			if (parent != null) {
				List<Integer> cousinLabels = this.getLabelSet(node == parent.getLeft() ? parent.getRight() : parent.getLeft());
				
				if (cousinLabels.contains(leftLabel)) {
					indicatorsInput.get().setValue(i, FlabelIndicator.RIGHT_BURST.val());
					return leftLabel;
				}
				
				if (cousinLabels.contains(rightLabel)) {
					indicatorsInput.get().setValue(i, FlabelIndicator.LEFT_BURST.val());
					return rightLabel;
				}
			}
			
			
			
			// Sample what type of burst this is
			double p1 = pInput.get().getValue(0); // Single mutant
			double p2 = pInput.get().getValue(1); // Double mutant
			if (p1 > 0 && p2 > 0) {
				p1 = 0.8; // Easier to initialise
				p2 = 0.2;
				//p2 = 0.5;
			}
			double u = Randomizer.nextDouble()*(p1 + p2);
			if (u < p1/2) {
				indicatorsInput.get().setValue(i, FlabelIndicator.LEFT_BURST.val());
				//node.setMetaData("INDICATOR", "LEFT_BURST");
				//node.setMetaData("LABEL", rightLabel);
				//Log.warning(node.getNr() + " has " + rightLabel);
				return rightLabel;
			}
			else if (u < p1) {
				indicatorsInput.get().setValue(i, FlabelIndicator.RIGHT_BURST.val());
				//node.setMetaData("INDICATOR", "leftLabel");
				//node.setMetaData("LABEL", rightLabel);
				//Log.warning(node.getNr() + " has " + leftLabel);
				return leftLabel;
			}
			else {
				indicatorsInput.get().setValue(i, FlabelIndicator.DOUBLE_BURST.val());
				nextLabelNr = nextLabelNr + 1;
				//node.setMetaData("INDICATOR", "DOUBLE_BURST");
				//node.setMetaData("LABEL", nextLabelNr);
				//Log.warning(node.getNr() + " has " + nextLabelNr);
				return nextLabelNr;
			}
			
		}
	
		
		
		
		
	}
	
	
	
	/**
	 * Log probability of a node, given its two children
	 * @param node
	 * @return
	 */
	public FlabelIndicator getIndicatorOfNode(Node node) {
		
		
		if (this.burstModelType == BurstModel.BINARY) {
			int val = indicatorsInput.get().getValue(node.getNr());
			return FlabelIndicator.get(val);
		}
		
		if (node.isLeaf()) {
			throw new IllegalArgumentException("Dev error: this method is for internal nodes only!");
		}
		int val = indicatorsInput.get().getValue(node.getNr()-this.nleaves);
		return FlabelIndicator.get(val);
	}
	
	
	
	public boolean isValid() {
		return update();
	}
	
	
	/**
	 * Note that internal node label ids are stochastic, but unique numbers
	 * @return
	 */
	public int[] getFlabels() {
		update();
		return this.nodeLabels;
	}

	
	/**
	 * Update internal node mapping from leave (including stubs), and indicators
	 */
	public boolean update() {
		
		//if (true) return true;
		
		// Todo efficiency - don't do this every time
		
		Tree tree = (Tree) treeInput.get();
		
		
		//return updateRecurse(tree.getRoot());
		int rootLabel = 0;
		nextLabelNr = 1;
		updateRecurse(tree.getRoot(), rootLabel);

		
		// Set leaf labels
		boolean valid = validateLeafLabels();
		
	
		
		// Reset
		if (!valid) {
			for (int leafNr = 0; leafNr < nleaves; leafNr++) {
				this.nodeLabels[leafNr] = getLeafLabel(tree.getNode(leafNr));
			}
		}
		
		return valid;
		
	}
	
	// Ensure that the new leaf labels form the same equivalence classes as the observed ones, even if the numbers are different
	private boolean validateLeafLabels() {
		
		if (simulatingInput.get()) return true;
		
		
		
		nextLabelNr = 0;
		
		// Equivalence classes. Each class contains all of the leaf indices with that true label value
		Map<Integer, List<Integer>> mapTrueLabels = new HashMap<>();
		for (int leafNr = 0; leafNr < nleaves; leafNr++) {
			int leafLabelTrue = getLeafLabel(treeInput.get().getNode(leafNr));
			if (mapTrueLabels.containsKey(leafLabelTrue)) {
				mapTrueLabels.get(leafLabelTrue).add(leafNr);
			}else {
				List<Integer> leafNrs = new ArrayList<>();
				leafNrs.add(leafNr);
				mapTrueLabels.put(leafLabelTrue, leafNrs);
			}
			
			nextLabelNr = Math.max(nextLabelNr, leafLabelTrue);
		}
		nextLabelNr++;
		
		
		// Same again but for the new labels
		Map<Integer, List<Integer>> mapNewLabels = new HashMap<>();
		for (int leafNr = 0; leafNr < nleaves; leafNr++) {
			int leafLabelNew = this.nodeLabels[leafNr];
			if (mapNewLabels.containsKey(leafLabelNew)) {
				mapNewLabels.get(leafLabelNew).add(leafNr);
			}else {
				List<Integer> leafNrs = new ArrayList<>();
				leafNrs.add(leafNr);
				mapNewLabels.put(leafLabelNew, leafNrs);
			}
		}
		
		
		
		// Now we need to confirm that the reassigned leaf labels produce the same classes, even if the labels are different
		Map<Integer, Integer> newToTrueMapping = new HashMap<>();
		for (int trueLabel : mapTrueLabels.keySet()) {
			
			List<Integer> leavesWithLabelTrue = mapTrueLabels.get(trueLabel);
			
			
			
			// Does there exist an equivalent list in the other mapping?
			for (int newLabel : mapNewLabels.keySet()) {
				List<Integer> leavesWithLabelNew = mapNewLabels.get(newLabel);
				
				
				boolean someAreTheSame = false;
				boolean someAreDifferent = false;
				for (int i : leavesWithLabelNew) {
					
					boolean contains = leavesWithLabelTrue.contains(i);
					if (contains) someAreTheSame = true;
					else someAreDifferent = true;
					
					if (someAreTheSame && someAreDifferent) break;
				}
				for (int i : leavesWithLabelTrue) {
					boolean contains = leavesWithLabelNew.contains(i);
					if (contains) someAreTheSame = true;
					else someAreDifferent = true;
					
					if (someAreTheSame && someAreDifferent) break;
				}
				
				// Not a valid assignment
				if (someAreTheSame && someAreDifferent) {
					
//					for (int t1 : mapTrueLabels.keySet()) {
//						List<Integer> list = mapTrueLabels.get(t1);
//						System.out.println("T" + t1 + ":" + Arrays.toString(list.toArray()));
//					}
//					
//					for (int t1 : mapNewLabels.keySet()) {
//						List<Integer> list = mapNewLabels.get(t1);
//						System.out.println("N" + t1 + ":" + Arrays.toString(list.toArray()));
//					}
//					
//					
//					
//					Log.warning("Invalid assignment for " + trueLabel + "/" + newLabel);
					return false;
				}
				
				if (someAreTheSame && !someAreDifferent) {
					
					// Found a match
					newToTrueMapping.put(newLabel, trueLabel);
					break;
					
					
				}
				
				
				
			}
			
			//if (foundMatch) break;
			
			
		}
		
	
		
		// Redo mapping
		int[] tmp = new int[this.nnodes];
		for (int nodeNr = 0; nodeNr < this.nnodes; nodeNr++) {
			int newLabel = this.nodeLabels[nodeNr];
			int trueLabel;
			if (newToTrueMapping.containsKey(newLabel)){
				trueLabel = newToTrueMapping.get(newLabel);
			}else {
				trueLabel = nextLabelNr;
				//Log.warning("next label " + nextLabelNr);
				nextLabelNr++;
			}
			
			tmp[nodeNr] = trueLabel;
			
		}
		
		for (int nodeNr = 0; nodeNr < this.nnodes; nodeNr++) {
			this.nodeLabels[nodeNr] = tmp[nodeNr];
		}
		
		//Log.warning("Good");
		//Log.warning("Valid assignment");
		return true;
	}
	
	
	// Pre order traversal
	private boolean updateRecurse(Node node, int labelThis) {
		
		
		// Any stubs on this branch?
		int nstubsOnBranch = stubsPerBranchInput.get() != null ? stubsPerBranchInput.get().getValue(node.getNr()) : stubsInput.get().getNStubsOnBranch(node.getNr());
		if (nstubsOnBranch > 0) {
			
			// New label
			labelThis = nextLabelNr;
			nextLabelNr++;
			
		}
		
		// Set label
		nodeLabels[node.getNr()] = labelThis;
		
		if (node.isLeaf()) {
			return true;
		}
		
		
		Node left = node.getLeft();
		Node right = node.getRight();
		
		
		// Get the indicator for this node 
		FlabelIndicator indicator = this.getIndicatorOfNode(node);
		switch (indicator) {
	
	
			case LEFT_BURST: {
				int labelLeft = nextLabelNr;
				nextLabelNr++;
				updateRecurse(left, labelLeft);
				
				int labelRight = labelThis;
				updateRecurse(right, labelRight);
				break;
			}
			
			case RIGHT_BURST: {
				
				int labelLeft = labelThis;
				updateRecurse(left, labelLeft);
				
				
				int labelRight = nextLabelNr;
				nextLabelNr++;
				updateRecurse(right, labelRight);
				break;
			}
			
			case DOUBLE_BURST: {
				
				int labelLeft = nextLabelNr;
				nextLabelNr++;
				updateRecurse(left, labelLeft);
				
				
				int labelRight = nextLabelNr;
				nextLabelNr++;
				updateRecurse(right, labelRight);
				break;
			}
			
			
			case NO_BURST: {
				updateRecurse(left, labelThis);
				updateRecurse(right, labelThis);
				break;
			}
			
		}
	
			
		
		
		return true;
	}
	
	
	private List<Integer> getLabelSet(Node node) {
		List<Integer> labels = new ArrayList<>();
		getLabelSet(node, labels);
		return labels;
	}
	
	
	
	private void getLabelSet(Node node, List<Integer> labels) {
		if (node.isLeaf()) {
			int label = getLeafLabel(node); 
			if (!labels.contains(label)) labels.add(label);
			return;
		}
		
		
		getLabelSet(node.getLeft(), labels);
		getLabelSet(node.getRight(), labels);
		
		
	}


	/**
	 *  Is this leaf's label represented at the root?
	 * @return
	 */
	private boolean isLeafLabelAtRoot(int nodeNr) {
		Tree tree = (Tree) treeInput.get();
		boolean val = isLeafLabelAtRootRecurse(tree.getNode(nodeNr));
		return val;
		
	}
	
	@Override
	protected void store() {
        super.store();
    }
	
	
	@Override
	protected void restore() {
        super.restore();
    }
	
	
	// Recurse upwards
	private boolean isLeafLabelAtRootRecurse(Node node) {
		
		
		if (node.isRoot()) return true;
		
		
		Node parent = node.getParent();
		FlabelIndicator parentIndicator = this.getIndicatorOfNode(parent);
		boolean thisIsLeftChild = node == parent.getLeft();
		
		// This child was NOT a single mutant
		if (parentIndicator == FlabelIndicator.DOUBLE_BURST || 
				(parentIndicator == FlabelIndicator.RIGHT_BURST && thisIsLeftChild) || 
				(parentIndicator == FlabelIndicator.LEFT_BURST && !thisIsLeftChild)) {
			return isLeafLabelAtRootRecurse(parent);
		}
		
		return false;
		
		
	}
	
	


	@Override
	public void init(PrintStream out) {
		
		
		out.print("nLabels\t");
		for (int i = 0; i < this.nleaves; i++) {
			String taxon = treeInput.get().getNode(i).getID();
			out.print("fleaf" + i + "\t");
		}
		for (int i = this.nleaves; i < this.nodeLabels.length; i++) {
			out.print("finternal" + i + "\t");
		}
		for (int i = 0; i < this.nleaves; i++) {
			String taxon = treeInput.get().getNode(i).getID();
			out.print("fleafAtRoot" + i + "\t");
		}
		
		
		
		out.print("nStubBursts\t");
		if (simulatingInput.get()) {
			for (int i = 0; i < this.nnodes; i++) {
				out.print("nburstsPerBranch." + (i+1) + "\t");
			}
		}
		
	}

	@Override
	public void log(long sample, PrintStream out) {
		
		int[] labels = this.getFlabels();
		
		// Number of labels
		List<Integer> labelsUnique = new ArrayList<>();
		for (int l : labels) {
			if (!labelsUnique.contains(l)) labelsUnique.add(l);
		}
		out.print(labelsUnique.size() + "\t");
		
		for (int i = 0; i < labels.length; i++) {
			out.print(labels[i] + "\t");
		}
		
		for (int i = 0; i < this.nleaves; i++) {
			boolean leafLabelAtRoot = this.isLeafLabelAtRoot(i);
			int val = leafLabelAtRoot ? 1 : 0;
			out.print(val + "\t");
		}
		
		if (simulatingInput.get()) {
			out.print(this.nBurstStubsTotal + "\t");
			for (int i = 0; i < this.nnodes; i++) {
				Node node= treeInput.get().getNode(i);
				out.print(nburstsPerBranch[node.getNr()] + "\t");
			}
		}else {
			
			int n = 0;
			for (int i = 0; i < this.treeInput.get().getNodeCount(); i++) {
				Node node = this.treeInput.get().getNode(i);
				n += getNumberOfBursts(node, true, false, false);
			}
			out.print(n + "\t");
			
		}
		
	}






	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
		
	}




	@Override
	public int getDimension() {
		
		if (this.burstModelType == BurstModel.BINARY) {
			return this.treeInput.get().getNodeCount();
		}else {
			return this.nodeLabels.length;
		}
		
		
	}




	@Override
	public double getArrayValue(int dim) {
		
		if (this.burstModelType == BurstModel.BINARY) {
			if (dim == this.treeInput.get().getRoot().getNr()) return 0;
			Node node = this.treeInput.get().getNode(dim);
			return this.getNumberOfBursts(node);
			
		}else {
			int[] labels = this.getFlabels();
			return labels[dim];
		}
		
		
	}
	
	
	

	
	

	@Override
    protected boolean requiresRecalculation() {

		
        if (InputUtil.isDirty(treeInput) || InputUtil.isDirty(leafLabelsInput) || InputUtil.isDirty(indicatorsInput)) {
       	 	return true;
        }
        
        if (stubsInput.get() != null && InputUtil.isDirty(stubsInput)) {
        	return true;
        }

        return false;
    }
	


	public int getNumberOfBursts(Node node) {
		return getNumberOfBursts(node, false, false, false);
	}
	
	
	public int getNumberOfSingleBursts(Node node) {
		return getNumberOfBursts(node, false, true, false);
	}
	
	public int getNumberOfDoubleBursts(Node node) {
		return getNumberOfBursts(node, false, false, true);
	}

	public int getNumberOfBursts(Node node, boolean stubsOnly, boolean singleOnly, boolean doubleOnly) {
		
		
		
		// Root does not have any changes from parent
		if (node.isRoot()) return 0;
		

		// Answer is precomputed
		if (simulatingInput.get()) {
			//if (true) return 0; 
			Log.warning("branch " + node.getNr() + " has " + nburstsPerBranch[node.getNr()] + " bursts");
			return nburstsPerBranch[node.getNr()];
		}
		
		
		
		
		
		
		//int[] labels = this.getFlabels();
		Node parent = node.getParent();
		
		//int n1 = node.getNr();
		//int n2 = parent.getNr();

		
		int nBursts = 0;

		
		if (!stubsOnly) {
			
			
			// Look at this node, not its parent
			if (this.burstModelType == BurstModel.BINARY) {
				FlabelIndicator indicator = getIndicatorOfNode(node);
				if (indicator == FlabelIndicator.DOUBLE_BURST) {
					nBursts++;
				}
			
			
			} else if (this.burstModelType != BurstModel.BINARY) {
				
			
				FlabelIndicator indicatorParent = getIndicatorOfNode(parent);
				boolean isLeftChild = node == parent.getLeft();
				
				
				if (singleOnly) {
					
					// Different label - yes there is a burst
					if ((indicatorParent == FlabelIndicator.RIGHT_BURST && !isLeftChild) || (indicatorParent == FlabelIndicator.LEFT_BURST && isLeftChild)) {
						nBursts++;
					}
					
				}
				
				else if (doubleOnly) {
					if (indicatorParent == FlabelIndicator.DOUBLE_BURST ) {
						nBursts++;
					}
				}
				
				else {
					
					// Different label - yes there is a burst
					if (indicatorParent == FlabelIndicator.DOUBLE_BURST 
							|| (indicatorParent == FlabelIndicator.RIGHT_BURST && !isLeftChild) 
							|| (indicatorParent == FlabelIndicator.LEFT_BURST && isLeftChild)) {
						nBursts++;
					}
					
				}
				
			}
		
			
			
			//Log.warning("bursts = " + nBursts + " for " + indicatorParent);
		
		}
		
		
		// How many stubs on this branch with changes
		for (int stubNr = 0; stubNr < this.getStubDimension(); stubNr++) {
			
			if (!this.stubs.includeStub(stubNr)) continue;
			if (this.stubs.getBranches().getValue(stubNr) != node.getNr()) continue;
			
			boolean stubIsAncestral = this.stubs.stubIsAncestral(stubNr);
			if (stubIsAncestral) {
				nBursts++;
			}
			
		}
		
		
		//Log.warning("bursts = " + nBursts + " for " + this.stubs.getNStubsOnBranch(node.getNr()));
		
		
		return nBursts;
	}
	


	public int getNStubs() {
		if (this.stubs == null) return 0;
		return this.stubs.getNStubs();
	}

	public int getStubDimension() {
		if (this.stubs == null) return 0;
		return this.stubs.getStubDimension();
	}


	public boolean includeStub(int stubNr) {
		if (this.stubs == null) return false;
		return this.stubs.includeStub(stubNr);
	}


	public int getIndicatorOfStub(int stubNr) {
		if (this.stubs == null) return 0;
		return this.stubs.getLabelIndicators().getValue(stubNr);
	}


	public boolean hasIndicators() {
		return this.stubs.hasIndicators();
	}



	public boolean isBurstNode(Node node) {
		FlabelIndicator indicator = getIndicatorOfNode(node);
		return indicator != FlabelIndicator.NO_BURST;
	}

	public int[] getLabels() {
		update();
		return nodeLabels;
	}


	

}




