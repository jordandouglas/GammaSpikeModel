package gammaspike.operator;

import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.operator.Uniform;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.StateNode;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import gammaspike.tree.Stubs;

@Description("Uniform operator on an internal node, aware of stubs")
public class StumpedTreeUniform extends Uniform {

	
	final public Input<Stubs> stubsInput = new Input<>("stubs", "stubs for the tree", Input.Validate.REQUIRED);
	
	
	
	@Override
    public double proposal() {
        final Tree tree = (Tree) InputUtil.get(treeInput, this);

        // randomly select internal node
        final int nodeCount = tree.getNodeCount();
        
        // Abort if no non-root internal nodes
        if (tree.getInternalNodeCount()==1) {
            return Double.NEGATIVE_INFINITY;
        }
        
        
        // Make sure that there is at least one non-fake and non-root internal node
        int leafNodeCount = tree.getLeafNodeCount();
        int fakeNodeCount = tree.getDirectAncestorNodeCount();
        if (fakeNodeCount == leafNodeCount-1 || (fakeNodeCount == leafNodeCount-2 && !tree.getRoot().isFake())) {
            return Double.NEGATIVE_INFINITY;
        }

        // Randomly select non-fake internal node
        Node node;
        do {
            final int nodeNr = nodeCount / 2 + 1 + Randomizer.nextInt(nodeCount / 2);
            node = tree.getNode(nodeNr);
        } while (node.isRoot() || node.isLeaf() || node.isFake());
        final double upper = node.getParent().getHeight();
        final double lower = Math.max(node.getLeft().getHeight(), node.getRight().getHeight());
        final double newValue = (Randomizer.nextDouble() * (upper - lower)) + lower;
       
        
        // Cache branch lengths before making proposal
        Stubs stubs = stubsInput.get();
        double[] cachedBranchLengths = stubs.prepareJacobian();
        
      
        
        // Make the change
        node.setHeight(newValue);
    
        
        // Jacobian. Relative stub heights stay the same but absolute heights change
        double logJacobian = stubs.getLogJacobian(cachedBranchLengths);
       

        return logJacobian;
    }
	 
	
	@Override
    public List<StateNode> listStateNodes() {
		Stubs stubs = stubsInput.get();
        final List<StateNode> list = super.listStateNodes();
       //if (stubAbsTimesConstantInput.get()) list.add(stubs.getStubHeights());
        list.add(stubs.getBranches());
        return list;
    }
	
}
