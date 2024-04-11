package gammaspike.clockmodel;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import gammaspike.flabel.Flabel;



@Description("Adds a burst of mutations after speciation events that result in a label change")
public class PunctuatedClockModel extends BranchRateModel.Base implements SpikeModel {
	
	final public Input<Tree> treeInput = new Input<>("tree", "the tree this relaxed clock is associated with.", Input.Validate.REQUIRED);
	final public Input<RealParameter> burstSizeInput = new Input<>("burstSize", "the additive clock rate after a label change.", Input.Validate.REQUIRED); 
	final public Input<Flabel> flabelsInput = new Input<>("flabel", "leaf and internal labels", Input.Validate.REQUIRED);
	final public Input<BooleanParameter> indicatorInput = new Input<>("indicator", "burst size is 0 of this is false", Input.Validate.OPTIONAL);
	
	
	int nnodes;
	double[] ratesArray;
	
	@Override
	public void initAndValidate() {
		
		this.nnodes = treeInput.get().getNodeCount();
        this.ratesArray = new double[this.nnodes];
	}

	
	public double getBurstSize() {
		if (indicatorInput.get() != null && !indicatorInput.get().getValue()) {
			return 0;
		}
		return burstSizeInput.get().getValue();
	}
	

	@Override
	public double getRateForBranch(Node node) {
		
		
		// Was there a change in function along this branch?
		Flabel flabels = flabelsInput.get();
		//flabels.update();
		
		int numberOfBursts = flabels.getNumberOfBursts(node);
		double burstRate = numberOfBursts * this.getBurstSize();
		double baseRate = meanRateInput.get().getArrayValue();
		double branchDistance = node.getLength()*baseRate + burstRate;
		
		//Log.warning((node.getNr()+1) + " has distance " + baseRate + "*" + node.getLength() + " + " + numberOfBursts + "*" + this.getBurstSize());
		
		
		if (node.isRoot()) return baseRate;
		
		// Effective rate takes into account burst and base rate
		double effectiveRate = branchDistance / node.getLength();
		
		//Log.warning(node.getID() + " has burst " + hasBurst + " rate=" + effectiveRate + " b = " + burstRate + " d= " + branchDistance);
		
		return effectiveRate;
	}

	
	@Override
    protected boolean requiresRecalculation() {

        if (InputUtil.isDirty(burstSizeInput) || InputUtil.isDirty(meanRateInput) || InputUtil.isDirty(flabelsInput)) {
       	 	return true;
        }
        
        
        if (indicatorInput.get() != null && InputUtil.isDirty(indicatorInput)) {
        	return true;
        }

        return false;
    }


	@Override
	public double[] getRatesArray() {
		for (int i = 0; i < nnodes; i ++) {
			Node node = this.treeInput.get().getNode(i);
    		ratesArray[i] = this.getRateForBranch(node);
    		//Log.warning("rate " + i + " = " + ratesArray[i]);
    	}
    	return ratesArray;
	}
	

}
