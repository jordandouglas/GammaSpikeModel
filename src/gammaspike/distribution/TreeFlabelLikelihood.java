package gammaspike.distribution;

import java.io.PrintStream;
import java.util.List;
import java.util.Random;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import gammaspike.flabel.Flabel;
import gammaspike.flabel.Flabel.BurstModel;
import gammaspike.flabel.Flabel.FlabelIndicator;


@Description("Probability of observing internal nod (including stub) labels, conidtional on leaf labels and tree")
public class TreeFlabelLikelihood extends Distribution implements Loggable {

	
	final public Input<TreeInterface> treeInput = new Input<>("tree", "tree without stubs", Input.Validate.REQUIRED);
	final public Input<Flabel> flabelsInput = new Input<>("flabel", "leaf and internal labels", Input.Validate.REQUIRED);
	final public Input<RealParameter> pInput = new Input<>("p", "the probability of a single-mutant", Input.Validate.REQUIRED);
	final public Input<Boolean> allowIllegalInitialisationInput = new Input<>("illegalInit", "allow the initial state to be illegal to help with initialisation", false);
	
	// TODO origin
	
	
	boolean initialising;
	
	@Override
	public void initAndValidate() {
		initialising = true;
		if (pInput.get().getDimension() != 1 && pInput.get().getDimension() != 3) {
			throw new IllegalArgumentException("Please ensure that p has 1 or 3 dimensions");
		}
		
		if (pInput.get().getDimension() == 3) {
			double p1 = pInput.get().getValue(0); // Single mutant
			double p2 = pInput.get().getValue(1); // Double mutant
			double p3 = pInput.get().getValue(2); // No mutant
			double sum = p1 + p2 + p3;
			if (Math.abs(1-sum) > 1e-6) {
				throw new IllegalArgumentException("Please ensure that p sums to 1, currently sums to " + sum);
			}
		}
		
	}
	
	private boolean initialisting() {
		return allowIllegalInitialisationInput.get() && initialising;
	}
	
	@Override
	public double calculateLogP() {
		
        logP = 0;
        
        
        // Load inputs
        Tree tree = (Tree) treeInput.get();
        
        int nEvents = tree.getInternalNodeCount(); 
        int nSinglesLeft = getNumBurstsOnInternalNodes(FlabelIndicator.LEFT_BURST);
        int nSinglesRight = getNumBurstsOnInternalNodes(FlabelIndicator.RIGHT_BURST);
        int ndoubleMutants = getNumBurstsOnInternalNodes(FlabelIndicator.DOUBLE_BURST);
        int nNoMutants = getNumBurstsOnInternalNodes(FlabelIndicator.NO_BURST);
        
        
        BurstModel burstModel = flabelsInput.get().getBurstModel();
        
        // Multinomial distribution. 
        // Do not include combinatorial component or it will double count equivalent states generated by the bit flip operator
        
        if (burstModel == BurstModel.BINARY) {
        	
        	
        	// Two outcomes - yes this branch has a spike or no it does not
        	double p = pInput.get().getValue();
        	logP = (ndoubleMutants)*Math.log(p) + (nNoMutants)*Math.log(1-p);
        	
        }
        
        
        
        else if (burstModel == BurstModel.THREE_STATE) {
        	
        	// Only three outcomes - left, right, double mutant
        	double p = pInput.get().getValue();
        	logP = (nSinglesLeft)*Math.log(p/2) + (nSinglesRight)*Math.log(p/2) + (ndoubleMutants)*Math.log(1-p);
        	
        }
        
        
        else if (burstModel == BurstModel.FOUR_STATE) {
        	
        	
        	// Ensure that the no-bursts only occur on valid trees
        	if (!flabelsInput.get().isValid()) {
        		//Log.warning("invalid - reject");
        		
        		logP = initialisting() ? -1e12 : Double.NEGATIVE_INFINITY;
        		initialising = false;
				return logP;
        	}
        	
        	
        	
        	
        	double p1 = pInput.get().getValue(0); // Single mutant
			double p2 = pInput.get().getValue(1); // Double mutant
			double p3 = pInput.get().getValue(2); // No mutant
			double sum = p1 + p2 + p3;
			if (Math.abs(1-sum) > 1e-6) {
				//Log.warning("reject1");
				logP = initialisting() ? -1e12 : Double.NEGATIVE_INFINITY;
				initialising = false;
				return logP;
			}
			
			if ((nSinglesLeft > 0 || nSinglesRight > 0) && p1 == 0) {
				//Log.warning("reject2");
				logP = initialisting() ? -1e12 : Double.NEGATIVE_INFINITY;
				initialising = false;
				return logP;
			}
			
			if (ndoubleMutants > 0 && p2 == 0) {
				//Log.warning("reject3");
				logP = initialisting() ? -1e12 : Double.NEGATIVE_INFINITY;
				initialising = false;
				return logP;
			}
			
			if (nNoMutants > 0 && p3 == 0) {
				
				//Log.warning("reject4");
				logP = initialisting() ? -1e12 : Double.NEGATIVE_INFINITY;
				initialising = false;
				return logP;
			}
			
			if (p1 > 0) {
				logP += (nSinglesLeft)*Math.log(p1/2) + (nSinglesRight)*Math.log(p1/2);
			}
			
			if (p2 > 0) {
				logP += (ndoubleMutants)*Math.log(p2);
			}
			
			if (p3 > 0){
				logP += (nNoMutants)*Math.log(p3);
			}
			
			
        	
        }

        initialising = false;
        return logP;
    }
	

	
	public int getNumBurstsOnInternalNodes(FlabelIndicator value) {
		
		int nBursts = 0;
		
        // Load inputs
        Flabel flabels = flabelsInput.get();
        Tree tree = (Tree) treeInput.get();
        int nnodes = tree.getNodeCount();
		int nleaves = tree.getLeafNodeCount();
        
        
        //flabels.update();
        
        // Iterate through internal nodes
        for (int nodeNr = nleaves; nodeNr < nnodes; nodeNr ++) {
        	Node node = tree.getNode(nodeNr);
        	if (node.isLeaf()) continue;
        	FlabelIndicator indicator = flabels.getIndicatorOfNode(node);
        	if (value == indicator) nBursts++;
        }
//        
//        
//        if (!includeStubs) return nSingles;
//       // nSingles = 0;
//        
//        if (!flabels.hasIndicators()) {
//        	return nSingles + flabels.getNStubs();
//        }
//        
//        // Iterate through stubs
//        for (int stubNr = 0; stubNr < flabels.getStubDimension(); stubNr++) {
//        	if (flabels.includeStub(stubNr)) {
//        		int indicator = flabels.getIndicatorOfStub(stubNr);
//        		if (indicator == +1 || indicator == -1) nSingles++;
//        	}
//        	
//        }
        
        return nBursts;
        
		
	}
	
	
	
	@Override
	public List<String> getArguments() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<String> getConditions() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void sample(State state, Random random) {
		// TODO Auto-generated method stub
		
	}
	
	
	
	@Override
    protected boolean requiresRecalculation() {
        return super.requiresRecalculation() || pInput.get().somethingIsDirty() ||  treeInput.get().somethingIsDirty() || InputUtil.isDirty(flabelsInput);
    }
	
	
	
	// Loggable
	@Override
	public void init(PrintStream out) {
		out.print("nSingleMutantsInternal\t");
	}
	
	@Override
	public void log(long sample, PrintStream out) {
		int nSingleMutantsInternal = getNumBurstsOnInternalNodes(FlabelIndicator.LEFT_BURST) + getNumBurstsOnInternalNodes(FlabelIndicator.RIGHT_BURST);
		out.print(nSingleMutantsInternal + "\t");
	}
	
	
	
	
	
	

}










