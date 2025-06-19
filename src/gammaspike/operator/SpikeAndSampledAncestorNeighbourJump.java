package gammaspike.operator;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;


@Description("Performs LeafToSampledAncestorJump on a leaf, and then adjusts the spike on the sibling branch")
public class SpikeAndSampledAncestorNeighbourJump extends TreeOperator {

	final public Input<RealParameter> spikesInput = new Input<>("spikes", "one spike size per branch.", Input.Validate.REQUIRED); 
	final public Input<Double> spikeMeanInput = new Input<>("spikeMean", "mean of the exponential distribuution for making spikes.", 0.001); 
	
	
    @Override
    public void initAndValidate() {
    }

    @Override
    public double proposal() {

        double newHeight, newRange, oldRange;

        Tree tree = treeInput.get();

        int leafNodeCount = tree.getLeafNodeCount();

        
        // What leaf
        Node leaf = tree.getNode(Randomizer.nextInt(leafNodeCount));
        Node parent = leaf.getParent();
        
        // Index of sibling
        if (parent.getChildCount() != 2) {
        	return Double.NEGATIVE_INFINITY; // Binary trees only
        }
        Node sibling = leaf == parent.getLeft() ? parent.getRight() : parent.getLeft();
        int siblingIndex = sibling.getNr();
        
        // Spike of the sibling
 		double spikeMean = spikeMeanInput.get();
 		double lambda = 1.0/spikeMean;
 		RealParameter spikes = spikesInput.get();
 		double sOldA = spikes.getValue(leaf.getNr());
 		double sOldB = spikes.getValue(siblingIndex);
 		double sNewA = 0;
 		double sNewB = 0;
 		
 		double pFwd, pRev;

        
        // We will go from (spikeA=0, lengthB=0, spikeB=0) to (spikeA>0, lengthB>0, spikeB>0), where A and B are siblings
        if (leaf.isDirectAncestor()) {
        	
        	if (sOldA != 0 || sOldB != 0) {
        		return Double.NEGATIVE_INFINITY; // Reject
        	}
        	
            oldRange = (double) 1;
            if (parent.isRoot()) {
                final double randomNumber = Randomizer.nextExponential(1);
                newHeight = parent.getHeight() + randomNumber;
                newRange = Math.exp(randomNumber);
            } else {
                newRange = parent.getParent().getHeight() - parent.getHeight();
                newHeight = parent.getHeight() + Randomizer.nextDouble() * newRange;
            }
            
            
            // Add a spike
            sNewA = Randomizer.nextExponential(lambda);
            sNewB = Randomizer.nextExponential(lambda);
			pRev = 0;
			pFwd = Math.log(lambda) - lambda*sNewA + Math.log(lambda) - lambda*sNewB; // Exponential density
            

        } 
        
        
     // We will go from (spikeA>0, lengthB>0) to (spikeA=0, lengthB=0), where A and B are siblings
        else {
        	
        	if (sOldA == 0 || sOldB == 0) {
        		return Double.NEGATIVE_INFINITY; // Reject
        	}
        	
        	
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
            
            
            // Delete the spike
            sNewA = 0;
            sNewB = 0;
			pFwd = 0;
			pRev = Math.log(lambda) - lambda*sOldB + Math.log(lambda) - lambda*sOldB; // Exponential density
			
        }
        
        // Update state
        spikes.setValue(siblingIndex, sNewA);
        spikes.setValue(leaf.getNr(), sNewB);
        parent.setHeight(newHeight);


        return Math.log(newRange/oldRange) + pRev - pFwd;
        
    }
    
    
    
    @Override
	public List<StateNode> listStateNodes() {
		final List<StateNode> list = new ArrayList<>();
		list.add(treeInput.get());
		list.add(spikesInput.get());
		return list;
   }

    
}
