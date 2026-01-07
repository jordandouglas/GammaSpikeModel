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


@Description("Performs LeafToSampledAncestorJump on a leaf, and then adjusts its spike")
public class SpikeAndSampledAncestorJump extends TreeOperator {

	final public Input<RealParameter> spikesInput = new Input<>("spikes", "one spike size per branch.", Input.Validate.REQUIRED); 
	final public Input<Double> spikeMeanInput = new Input<>("spikeMean", "mean of the exponential distribution for making spikes.", 0.001);
	
	
    @Override
    public void initAndValidate() {
    }

    @Override
    public double proposal() {

        double newHeight, newRange, oldRange;

        Tree tree = treeInput.get();

        int leafNodeCount = tree.getLeafNodeCount();


        // Randomly select a leaf node
        Node leaf = tree.getNode(Randomizer.nextInt(leafNodeCount));
        Node parent = leaf.getParent();
        
        
        // Spike of the sampled ancestor
 		double spikeMean = spikeMeanInput.get();
 		double lambda = 1.0/spikeMean;
 		RealParameter spikes = spikesInput.get();
 		double sOld = spikes.getValue(leaf.getNr());
 		double sNew = 0;
 		double pFwd, pRev;

        
        // We will go from (spike=0, length=0) to (spike>0, length>0)
        if (leaf.isDirectAncestor()) {
        	
        	if (sOld != 0) {
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
            sNew = Randomizer.nextExponential(lambda);
			pRev = 0;
			pFwd = Math.log(lambda) - lambda*sNew; // Exponential density
            

        } 
        
        
        // We will go from (spike>0, length>0) to (spike=0, length=0)
        else {
        	
        	if (sOld == 0) {
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
            sNew = 0;
			pFwd = 0;
			pRev = Math.log(lambda) - lambda*sOld; // Exponential density
			
        }
        
        // Update state
        spikes.setValue(leaf.getNr(), sNew);
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
