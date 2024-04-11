package gammaspike.operator;

import java.util.ArrayList;
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
import gammaspike.tree.Stubs.Stub;

// BROKEN operator do not use
@Description("Moves a stub from one branch to another at the same height")
public class StubBranchOperator extends Operator {
	
	
	final public Input<Stubs> stubsInput = new Input<>("stubs", "stubs for the tree", Input.Validate.REQUIRED);
	final public Input<TreeInterface> treeInput = new Input<>("tree", "the tree", Input.Validate.REQUIRED);

	
	
	@Override
	public void initAndValidate() {
		
		
	}

	
	
	
	@Override
	public double proposal() {
		
		// Does not work! Hastings ratio is wrong
		if (true) return Double.NEGATIVE_INFINITY;
		
		// Get parameters
		Stubs stubs = stubsInput.get();
		RealParameter times = stubs.getStubHeights();
		IntegerParameter branches = stubs.getBranches();
		IntegerParameter labelIndicators = stubs.getLabelIndicators();
		Tree tree = (Tree)treeInput.get();
		int nnodes = treeInput.get() != null ? tree.getNodeCount()-1 : 1; // Origin not supported

		if (stubs.getNStubs() == 0) {
			return Double.NEGATIVE_INFINITY;
		}
		
		// Pick a stub
		Stub stubToMove = stubs.getRandomStub();
		int nodeFromNr = stubToMove.getBranchNr();
		double branchLengthOriginal = tree.getNode(nodeFromNr).getLength();
		
		// Pick a destination branch
		int nodeToNr = Randomizer.nextInt(nnodes);
		while (nodeToNr == nodeFromNr) {
			nodeToNr = Randomizer.nextInt(nnodes);
		}
		double branchLengthNew = tree.getNode(nodeToNr).getLength();
		
		
		
		// Move it
		branches.setValue(stubToMove.getIndex(), nodeToNr);
		
		//Log.warning("Moving from " + nodeFromNr + " to " + nodeToNr + " with lengths " + branchLengthOriginal+ " and " + branchLengthNew);
		
		// Hastings ratio is length ratio
		return Math.log(branchLengthNew/branchLengthOriginal);
		
//		
//		// Does not work! Hastings ratio is wrong
//		//if (true) return Double.NEGATIVE_INFINITY;
//		
//		
//		// Get parameters
//		Stubs stubs = stubsInput.get();
//		RealParameter times = stubs.getStubHeights();
//		IntegerParameter branches = stubs.getBranches();
//		IntegerParameter labelIndicators = stubs.getLabelIndicators();
//		Tree tree = (Tree)treeInput.get();
//		int nnodes = treeInput.get() != null ? tree.getNodeCount()-1 : 1; // Origin not supported
//	
//		
//		if (stubs.getNStubs() == 0) {
//			return Double.NEGATIVE_INFINITY;
//		}
//		
//		// Pick a stub
//		Stub stubToMove = stubs.getRandomStub();
//		int nodeFromNr = stubToMove.getBranchNr();
//		double stubHeight = stubToMove.getAbsoluteHeight();
//		
//		
//		
//		// Get all branches that cross this height
//		List<Integer> branchesAtRightHeight = new ArrayList<>();
//		for (int branchNr = 0; branchNr < tree.getNodeCount(); branchNr ++) {
//			Node node = tree.getNode(branchNr);
//			if (branchNr == nodeFromNr) continue;
//			if (node.isRoot()) continue;
//
//			double t0 = node.getHeight();
//			double t1 = node.getParent().getHeight();
//			
//			if (t0 <= stubHeight && t1 > stubHeight) {
//				branchesAtRightHeight.add(branchNr);
//			}
//		}
//		if (branchesAtRightHeight.isEmpty()) return Double.NEGATIVE_INFINITY;
//		
//		
//		// Sample a destination branch
//		int nodeToNr = Randomizer.nextInt(branchesAtRightHeight.size());
//		int n2 = stubs.getNStubsOnBranch(nodeToNr);
//		double newRelativeHeight = (stubHeight - tree.getNode(nodeToNr).getHeight()) / tree.getNode(nodeToNr).getLength();
//		
//		
//		
//		// Move the stub to the other branch. Its absolute height will be the same but relative height will change
//		branches.setValue(stubToMove.getIndex(), nodeToNr);
//		times.setValue(stubToMove.getIndex(), newRelativeHeight);
//		
//		
//		// How many branches have stubs now?
//		//List<Integer> branchesWithStubsAfter = stubs.getBranchesWithStubs();
//
//		
//		return 0;
		
		
	}
	
	
	
	
	@Override
    public List<StateNode> listStateNodes() {
		Stubs stubs = stubsInput.get();
        final List<StateNode> list = super.listStateNodes();
       // list.add(stubs.getStubHeights());
        list.add(stubs.getBranches());
        return list;
    }

}





