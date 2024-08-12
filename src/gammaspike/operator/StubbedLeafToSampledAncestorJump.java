package gammaspike.operator;

import beast.base.core.Description;
import beast.base.core.Input;
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
