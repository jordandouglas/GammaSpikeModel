package gammaspike.operator;

import java.util.List;

import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.inference.Operator;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import gammaspike.distribution.StumpedTreePriorOneLineage;
import gammaspike.tree.Stubs;

public class StubGibbsOperatorSingleLineage extends Operator {
	
	
	final public Input<StumpedTreePriorOneLineage> priorInput = new Input<>("prior", "stubs for the tree", Input.Validate.REQUIRED);
	
	Stubs stubs;

	@Override
	public void initAndValidate() {
		this.stubs = priorInput.get().stubsInput.get();
		if (!this.stubs.getReversibleJump()) {
			throw new IllegalArgumentException("Cannot use this operator unless reversible jump is enabled");
		}
	}

	
	@Override
	public double proposal() {
		
		
		
		// Sample a branch
		StumpedTreePriorOneLineage prior = priorInput.get();
		
		
		// What is the prior of this branch before the proposal?
		double logPBranchBefore = prior.calculateLogP();
		double totalTime = prior.getTotalTime();
		
		
		// Store
		stubs.storeDimensions();
		
		
		// Sample number of stubs on this branch
		double g = prior.getExpectedNumStubs();
		int oldM = stubs.getNStubs();
		int newM = (int) Randomizer.nextPoisson(g);
		
		
		boolean print = false; //oldM != newM;
		
		// Change dimension
		RealParameter times = stubs.getStubHeights();
		IntegerParameter branches = stubs.getBranches();
		times.setDimension(newM+1);
		branches.setDimension(newM+1);
		
		
		
		if (print) Log.warning("Proposed M from " + oldM + " to " + stubs.getNStubs());
		
		// Resample all values
		for (int i = 1; i <= newM; i ++) {
			
			
			double oldTime = times.getValue(i);
			double time = prior.sampleTime();
			
			// Normalise between 0 and 1
			time = totalTime - time;
			time = time / totalTime;
			times.setValue(i, time);
			
			if (print) Log.warning("Proposed time from " + oldTime + " to " + time);
			
			
			
		}
		
		
		
		// What is the prior of this branch after the proposal?
		double logPBranchAfter = prior.calculateLogP();
		
		
		if (print) Log.warning("pbefore " + logPBranchBefore + " pafter " + logPBranchAfter);
		
		return logPBranchBefore - logPBranchAfter;
		
	}
	
	

	@Override
	public void accept() {
		stubs.acceptDimensions();
		super.accept();
	}
	
	@Override
	public void reject(final int reason) {
		super.reject(reason);
	}
	
	@Override
    public List<StateNode> listStateNodes() {
        final List<StateNode> list = super.listStateNodes();
        list.add(stubs.getStubHeights());
        list.add(stubs.getBranches());
        return list;
    }
	
	

}
