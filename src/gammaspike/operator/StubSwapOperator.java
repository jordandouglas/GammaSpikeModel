package gammaspike.operator;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.Operator;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.Parameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;

@Description("Swap operator but it rejects proposal if there is just 1 element")
public class StubSwapOperator extends Operator {
	
	
    final public Input<IntegerParameter> intparameterInput = new Input<>("parameter", "an integer parameter to swap individual values for", Validate.REQUIRED);

    IntegerParameter parameter;
    
    @Override
    public void initAndValidate() {
        parameter = intparameterInput.get();
    }
    
    
    @Override
    public double proposal() {
    	
    	if (parameter.getDimension() < 3) {
    		return Double.NEGATIVE_INFINITY;
    	}
    	
    	
    	// Skip dummy at index 0
    	List<Integer> allIndices = new ArrayList<>();
    	for (int i = 1; i < parameter.getDimension(); i++) {
    		allIndices.add(i);
        }
    	
    	
        int left, right;
        for (int i = 0; i < 1; i++) {
            left = allIndices.remove(Randomizer.nextInt(allIndices.size()));
            right = allIndices.remove(Randomizer.nextInt(allIndices.size()));
            parameter.swap(left, right);
        }

        return 0.0;
    }





}
