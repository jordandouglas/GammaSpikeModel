package gammaspike.operator;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.operator.TreeStretchOperator;
import gammaspike.tree.Stubs;



@Description("Tree stretchflex operator that accounts for stub placement")
public class StubTreeStretchOperator extends TreeStretchOperator {
	

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
