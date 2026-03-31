package gammaspike.operator;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.spec.evolution.operator.UpDownOperator;
import gammaspike.tree.Stubs;


// Cross between scale and updown
@Description("Scales nodes and stubs in a tree")
public class StumpedTreeScaler extends UpDownOperator {


	final public Input<Stubs> stubsInput = new Input<>("stubs", "stubs for the tree", Input.Validate.OPTIONAL);
	

    @Override
    public void initAndValidate() {
    	super.initAndValidate();
    }


    /**
     * Same proposal as UpDown, but we must also account for stubs HR if they are being estimated directly
     */
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

