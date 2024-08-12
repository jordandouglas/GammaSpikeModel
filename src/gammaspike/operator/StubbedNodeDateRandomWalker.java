package gammaspike.operator;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.operator.TipDatesRandomWalker;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;
import gammaspike.tree.Stubs;
import sa.evolution.operators.SampledNodeDateRandomWalker;
import sa.evolution.tree.SamplingDate;
import sa.evolution.tree.TreeWOffset;



// Same as SampledNodeDateRandomWalker (as package) but accounts for stubs

@Description("Randomly select a sampled node and shifts the date of the node within a given window")
public class StubbedNodeDateRandomWalker extends SampledNodeDateRandomWalker {
	
	
	final public Input<Stubs> stubsInput = new Input<>("stubs", "stubs for the tree", Input.Validate.REQUIRED);

	
	@Override
    public double proposal() {
		
		// Cache branch lengths before making proposal
		Stubs stubs = stubsInput.get();
        double[] cachedBranchLengths = stubs.prepareJacobian();
		
		double HR = super.proposal();
		
		// Jacobian. Relative stub heights stay the same but absolute heights change
        double logJacobian = stubs.getLogJacobian(cachedBranchLengths);
        
        return HR + logJacobian;
		
	}


}
