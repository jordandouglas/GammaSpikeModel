package gammaspike.operator;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;
import gammaspike.tree.Stubs;
import sa.evolution.operators.LeafToSampledAncestorJump;


@Description("Same as LeafToSampledAncestorJump but with stub awareness")
public class StubbedLeafToSampledAncestorJump extends LeafToSampledAncestorJump  {
    
	final public Input<Stubs> stubsInput = new Input<>("stubs", "stubs for the tree", Input.Validate.REQUIRED);
	
	
	boolean rootMove = false;
	
	@Override
    public double proposal() {
		
		// Cache branch lengths before making proposal
		Stubs stubs = stubsInput.get();
        double[] cachedBranchLengths = stubs.prepareJacobian();
		
		double HR = proposaInnerl(); //super.proposal();
		
		// Jacobian. Relative stub heights stay the same but absolute heights change
        double logJacobian = stubs.getLogJacobian(cachedBranchLengths);
        
        return HR + logJacobian;
		
	}
	
	
    public double proposaInnerl() {

        double newHeight, newRange, oldRange;
        int categoryCount = 1;
        if (categoriesInput.get() != null) {

            categoryCount = categoriesInput.get().getUpper() - categoriesInput.get().getLower() +1;
        }

        Tree tree = treeInput.get();

        int leafNodeCount = tree.getLeafNodeCount();

        Node leaf = tree.getNode(Randomizer.nextInt(leafNodeCount));
        Node parent = leaf.getParent();
        

        if (leaf.isDirectAncestor()) {
            oldRange = (double) 1;
            if (parent.isRoot()) {
            	
                final double randomNumber = Randomizer.nextExponential(10.0/parent.getHeight()); // Use a more informed jump size (not just 1)
                newHeight = parent.getHeight() + randomNumber;
                newRange = Math.exp(randomNumber);
                
            } else {
            	
                newRange = parent.getParent().getHeight() - parent.getHeight();
                newHeight = parent.getHeight() + Randomizer.nextDouble() * newRange;
            }

            if (categoriesInput.get() != null) {
                int index = leaf.getNr();
                int newValue = Randomizer.nextInt(categoryCount) + categoriesInput.get().getLower(); // from 0 to n-1, n must > 0,
                categoriesInput.get().setValue(index, newValue);
            }
        } else {
        	
        	
        	
            newRange = (double) 1;
            //make sure that the branch where a new sampled node to appear is not above that sampled node
            if (getOtherChild(parent, leaf).getHeight() >= leaf.getHeight())  {
                return Double.NEGATIVE_INFINITY;
            }
            if (parent.isRoot()) {
                oldRange = Math.exp(parent.getHeight() - leaf.getHeight());
            } else {
                oldRange = parent.getParent().getHeight() - leaf.getHeight();
            }
            newHeight = leaf.getHeight();
            if  (categoriesInput.get() != null) {
                int index = leaf.getNr();
                categoriesInput.get().setValue(index, -1);
            }
        }
        parent.setHeight(newHeight);
        
        

        //make sure that either there are no direct ancestors or r<1
        if ((rInput.get() != null) && (tree.getDirectAncestorNodeCount() > 0 && rInput.get().getValue() == 1))  {
            return Double.NEGATIVE_INFINITY;
        }

        return Math.log(newRange/oldRange);
    }
    

}
