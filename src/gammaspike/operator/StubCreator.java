package gammaspike.operator;

import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.Operator;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import gammaspike.tree.Stubs;


@Description("Create or destroy a stub")
public class StubCreator extends Operator {
	
	
	final public Input<Stubs> stubsInput = new Input<>("stubs", "stubs for the tree", Input.Validate.REQUIRED);
	final public Input<RealParameter> totalTimeInput = new Input<>("totalTime", "totalTime", Input.Validate.OPTIONAL);
	final public Input<TreeInterface> treeInput = new Input<>("tree", "the tree", Input.Validate.XOR, totalTimeInput);
	
	
	final double DIRECT_ANCESTOR_THRESHOLD = 1e-12;
	Stubs stubs;
	int nextNr = 0;

	@Override
	public void initAndValidate() {
		this.stubs = stubsInput.get();
		if (!this.stubs.getReversibleJump()) {
			throw new IllegalArgumentException("Cannot use this operator unless reversible jump is enabled");
		}
	}

	@Override
	public double proposal() {
		

		RealParameter times = stubs.getStubHeights();
		IntegerParameter branches = stubs.getBranches();
		IntegerParameter labelIndicators = stubs.getLabelIndicators();
		
		//Log.warning("upper = " + branches.getUpper());
		
		Tree tree = (Tree) treeInput.get();
		int nnodes = treeInput.get() != null ? tree.getNodeCount()-1 : 1; // Origin not supported
	
		
		
		// Sample a branch first
		int nodeNr = Randomizer.nextInt(nnodes);
		Node node = tree.getNode(nodeNr);
		while (node.isLeaf() && (node.getParent().isDirectAncestor() || node.getLength() <= DIRECT_ANCESTOR_THRESHOLD)) {
			nodeNr = Randomizer.nextInt(nnodes);
			node = tree.getNode(nodeNr);
		}
		
		// Sample indicator
		int indicatorNr = Randomizer.nextInt(3) - 1;
		
		//nextNr ++;
		//if (nextNr > 1) nextNr = 0;
		//Log.warning("indicatorNr " + indicatorNr);
		
		// Test
		//nodeNr = 3;
		//nnodes = 1;
		
		double dt = getDt(nodeNr);
		int oldM = stubs.getNStubs();
		int oldMBranchAllIndicators = stubs.getNStubsOnBranch(nodeNr);
		int oldMBranch = oldMBranchAllIndicators;
		if (labelIndicators != null) {
			oldMBranch = stubs.getNStubsOnBranchWithIndicator(nodeNr, indicatorNr);
			//Log.warning(oldMBranch + " on branch with indicator " + indicatorNr + " out of " + stubs.getNStubsOnBranch(nodeNr));
		}
		
		
		//double h = tree.getNode(nodeNr).getHeight();
		
		//Log.warning("on branch " + nodeNr + " dt=" + dt);
		
		
		// The effective dimension is this number - 1 (the first element is a dummy placeholder)
		int dimension = stubs.getStubDimension();
		
		// Add or remove a dimension?
		boolean create = Randomizer.nextBoolean();
		if (create) {
			
			
			// Make a new stub on this branch
			int newMBranch = oldMBranch + 1;
			int newM = oldM + 1;
			
			//Log.warning("Going from " + oldMBranch + " to " + newMBranch + " for node " + nodeNr +  " out of " + nnodes + " with i=" + indicatorNr);
			
			
			//Log.warning("StubCreator add");
			stubs.storeDimensions();
			
			//Log.warning("CREATE: from " + stubs.getNStubs());
			
			
			// Create new element
			int newDimension = dimension+1;
			times.setDimension(newDimension);
			branches.setDimension(newDimension);
			if (labelIndicators != null) labelIndicators.setDimension(newDimension);
			

			// Add new value to a random position in the list, and shift everything after it downstream
			int indexToInsertAt = dimension == 1 ? 1 : Randomizer.nextInt(dimension-1)+1; // Insert at 1,2,...,d. Do not disrupt the one at 0
			for (int oldIndex = dimension-1; oldIndex >= indexToInsertAt; oldIndex--) {
				int newIndex = oldIndex+1;
				times.setValue(newIndex, times.getValue(oldIndex));
				branches.setValue(newIndex, branches.getNativeValue(oldIndex));
				if (labelIndicators != null) {
					labelIndicators.setValue(newIndex, labelIndicators.getNativeValue(oldIndex));
				}
			}
			
			
			double u = Randomizer.nextDouble(); 
			
			
			// Set new values
			times.setValue(indexToInsertAt, u); // Relative time along branch
			branches.setValue(indexToInsertAt, nodeNr);
			if (labelIndicators != null) {
				labelIndicators.setValue(indexToInsertAt, indicatorNr);
			}

			double logHR = Math.log(dt); 
			if (labelIndicators != null) {
				//logHR += Math.log(3); // Three possible values for new label: -1, 0, +1
				logHR += Math.log(oldMBranchAllIndicators+1) - Math.log(newMBranch) - Math.log(3);
			}
			
			
			return logHR;
			
			
			
			
		}else {
			
			
			// Delete a stub from this branch and with this indicator
			int newMBranch = oldMBranch - 1;
			int newM = oldM - 1;
			if (newMBranch < 0) return Double.NEGATIVE_INFINITY;
			
			
			//Log.warning("Going from " + oldMBranch + " to " + newMBranch + " for node " + nodeNr + " with i=" + indicatorNr);
			
			stubs.storeDimensions();
			


			int newDimension = dimension-1;
			
			 // Destroy one from index 1,2,...,d such that it is on the right branch
			int indexToDestroyOnBranch = Randomizer.nextInt(oldMBranch);
			
			
			// Find the corresponding stub index on this branch
			int indexToDestroy = -1;
			int stubNrOnBranch = 0;
			for (int i = 0; i < dimension; i ++) {
				if (!stubs.includeStub(i)) continue;
				int b = tree != null ? stubs.getBranches().getValue(i) : nodeNr;
				boolean indicatorMatch = labelIndicators == null ? true : stubs.getLabelIndicatorOfStub(i) == indicatorNr;
				if (indicatorMatch && b == nodeNr) {
					if (stubNrOnBranch == indexToDestroyOnBranch) {
						indexToDestroy = i;
						break;
					}
					stubNrOnBranch++;
				}
			}
			
			
			if (indexToDestroy == -1) {
				Log.warning("Unexpected error " + indexToDestroy + " _ " + indexToDestroyOnBranch + " / " + oldMBranch);
			}
			
			// Shift everything after the deletion point down by 1
			for (int oldIndex = indexToDestroy+1; oldIndex < dimension; oldIndex++) {
				
				int newIndex = oldIndex-1;
				times.setValue(newIndex, times.getValue(oldIndex));
				branches.setValue(newIndex, branches.getNativeValue(oldIndex));
				if (labelIndicators != null) {
					labelIndicators.setValue(newIndex, labelIndicators.getNativeValue(oldIndex));
				}
			}
			
			// Remove the last entry
			times.setDimension(newDimension);
			branches.setDimension(newDimension);
			if (labelIndicators != null) {
				labelIndicators.setDimension(newDimension);
			}
			
			// Changing dimension resets the dirtiness so now we set it to dirty again
			for (int i = 0; i < newDimension; i++) {
				times.setValue(i, times.getValue(i));
				branches.setValue(i, branches.getNativeValue(i));
				if (labelIndicators != null) {
					labelIndicators.setValue(i, labelIndicators.getNativeValue(i));
				}
			}
			//Log.warning("destroy dt = " + dt);
			
			
			double logHR = -Math.log(dt); 
			if (labelIndicators != null) {
				//logHR -= Math.log(3); // Three possible values for new label: -1, 0, +1
				logHR += Math.log(oldMBranch) + Math.log(3) - Math.log(oldMBranchAllIndicators);
			}
			

			return logHR;
			
			
		}
		

	

		
		
	}

	
	/**
	 * The length of the lineage. Total time if just 1 lineage, selected branch length if tree, 1 otherwise
	 * @param branchNr
	 * @return
	 */
	private double getDt(int branchNr) {
		
		
		if (totalTimeInput.get() != null) {
			return totalTimeInput.get().getValue();
		}
		
		
		else if (treeInput.get() != null) {
			Tree tree = (Tree) treeInput.get();
			Node node = tree.getNode(branchNr);
			if (node.isRoot()) return 0; // TODO: ORIGIN NOT SUPPORTED
			return node.getParent().getHeight() - node.getHeight();
		}
		
		return 1;
		
		
		
	}
	
	
	@Override
	public void accept() {
		stubs.acceptDimensions();
		super.accept();
	}
	
	@Override
	public void reject(final int reason) {
		super.reject(reason);
	}
	
	@Override
    public List<StateNode> listStateNodes() {
        final List<StateNode> list = super.listStateNodes();
        list.add(stubs.getStubHeights());
        list.add(stubs.getBranches());
        if (stubs.getLabelIndicators() != null) list.add(stubs.getLabelIndicators());
        return list;
    }
	
	
	
}
