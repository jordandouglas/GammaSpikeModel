package gammaspike.operator;

import java.util.List;

import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Operator;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import gammaspike.distribution.StumpedTreePrior;
import gammaspike.tree.Stubs;

public class StubGibbsOperator extends Operator {
	
	
	final public Input<StumpedTreePrior> priorInput = new Input<>("prior", "stubs for the tree", Input.Validate.REQUIRED);
	final public Input<Boolean> incrementInput = new Input<>("increment", "change number of stubs by +/- 1 only", Input.Validate.REQUIRED);
	
	Stubs stubs;

	@Override
	public void initAndValidate() {
		this.stubs = priorInput.get().stubsInput.get();
		if (!this.stubs.getReversibleJump()) {
			throw new IllegalArgumentException("Cannot use this operator unless reversible jump is enabled");
		}
	}

	
	@Override
	public double proposal() {
		
		StumpedTreePrior prior = priorInput.get();
		

		double totalTime = prior.getTotalTime();
		
		// Sample a branch
		Tree tree = prior.getTree();
		int nodeNr = Randomizer.nextInt(tree.getNodeCount()-1); // Excluding root
		double nodeHeight = tree.getNode(nodeNr).getHeight();
		double nodeParentHeight = tree.getNode(nodeNr).getParent().getHeight();
		double nodeTime = totalTime - nodeHeight;
		double nodeParentTime = totalTime - nodeParentHeight;
		
		
		
		// What is the prior of this branch before the proposal?
		double logPBranchBefore = prior.getLogPForBranch(nodeNr);
		
		
		
		// Sample number of stubs on this branch
		double g = prior.getGForInterval(nodeHeight, nodeParentHeight);
		int oldM = stubs.getNStubs();
		int oldMBranch = stubs.getNStubsOnBranch(nodeNr);
		int newMBranch = (int) Randomizer.nextPoisson(g);
		
		if (incrementInput.get()) {
			newMBranch = oldMBranch + (Randomizer.nextBoolean() ? -1 : +1);
			if (newMBranch < 0) return Double.NEGATIVE_INFINITY;
		}
		int newM = oldM + newMBranch - oldMBranch;
		
		// Nothing changes
		if (oldMBranch == 0 && newMBranch == 0) {
			return 0;
		}
		
		// Store
		stubs.storeDimensions();
		
		
		boolean print = false; //oldM != newM;
		
		// Change dimension
		RealParameter heights = stubs.getStubHeights();
		IntegerParameter branches = stubs.getBranches();


		
		if (print) Log.warning("Proposed M from " + oldMBranch + " to " + newMBranch + " for branch " + nodeNr);
		
		
		//if (print) Log.warning("newM = " + newM + " oldM = " + oldM + " newMBranch=" + newMBranch + " oldMBranch=" + oldMBranch);
		
		
		// For easier bookkeeping, shift all other branch entries down to the start of the list
		int[] newBranchNrs = new int[newM+1];
		double[] newHeights = new double[newM+1];
		int newIndex = 0;
		for (int i = 0; i <= oldM; i ++) {
			int branchNr = branches.getValue(i);
			double height = heights.getValue(i);
			if (i == 0 || branchNr != nodeNr) {
				newBranchNrs[newIndex] = branchNr;
				newHeights[newIndex] = height;
				newIndex++;
			}
		}
		
		
		//if (print) Log.warning("Relocated " + newIndex + " / " + newM);
		
		
		// Debugging
		int expectedNumberOfProposals = newMBranch;
		int numproposals = 0;
		
		// Now, resample the ones at the end of the list (if there are any)
		for (int i = newIndex; i <= newM; i ++) {
			
			double time = prior.sampleTime(nodeParentTime, nodeTime);
			
			// Convert from time to height
			double height = totalTime - time;
			
			// Normalise height between 0 and 1
			height = height / totalTime;
			
			//if (print) Log.warning("Proposed height " + height);
			
			newBranchNrs[newIndex] = nodeNr;
			newHeights[newIndex] = height;
			newIndex++;
			
			numproposals++;
			
		}
		
		if (newIndex != newHeights.length) {
			Log.warning("Proposed M from " + oldMBranch + " to " + newMBranch + " for branch " + nodeNr);
			Log.warning("Only populated " + newIndex + " / " + newHeights.length + " newM=" + newM);
			System.exit(1);
		}
		
		
		// Sanity checks
		if (expectedNumberOfProposals != numproposals) {
			Log.warning("Proposed M from " + oldMBranch + " to " + newMBranch + " for branch " + nodeNr);
			Log.warning("Expected " + expectedNumberOfProposals + " proposals but made " + numproposals);
			System.exit(1);
		}
		
		
		// Copy over to parameters array
		heights.setDimension(newM+1);
		branches.setDimension(newM+1);
		for (int i = 0; i <= newM; i ++) {
			heights.setValue(i, newHeights[i]);
			branches.setValue(i, newBranchNrs[i]);
		}
		
		// What is the prior of this branch after the proposal?
		double logPBranchAfter = prior.getLogPForBranch(nodeNr);
		
		
		
		
		
		
		// Adjust hastings ratios by removing poisson term
		double poissonAdjust = 0;
		if (incrementInput.get()) {
			
			double mAdjust = newMBranch > oldMBranch ? Math.log(newMBranch) : -Math.log(oldMBranch);
			mAdjust = mAdjust * -1;
			double g2 = (newMBranch - oldMBranch) * Math.log(g);
			poissonAdjust = mAdjust + g2;
			
		}
		
		if (print) Log.warning("pbefore " + logPBranchBefore + " pafter " + logPBranchAfter + " adjust " + poissonAdjust);
		if (print) Log.warning("HR = " + (logPBranchBefore - logPBranchAfter + poissonAdjust));
		
		return logPBranchBefore - logPBranchAfter + poissonAdjust;
		
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
        return list;
    }
	
	

}
