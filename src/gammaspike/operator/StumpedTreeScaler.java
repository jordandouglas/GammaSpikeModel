package gammaspike.operator;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Operator;
import beast.base.inference.StateNode;
import beast.base.inference.operator.kernel.KernelDistribution;
import beast.base.inference.parameter.Parameter;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import gammaspike.tree.Stubs;


// Cross between scale and updown
@Description("Scales nodes and stubs in a tree")
public class StumpedTreeScaler extends Operator {


	final public Input<Tree> treeInput = new Input<>("tree", "all tree divergence times are scaled", Input.Validate.REQUIRED);
	final public Input<Stubs> stubsInput = new Input<>("stubs", "stubs for the tree", Input.Validate.REQUIRED);
	final public Input<List<StateNode>> upInput = new Input<>("up", "zero or more items to scale in the same direction as the tree", new ArrayList<>());
	final public Input<List<StateNode>> downInput = new Input<>("down", "zero or more items to scale in the reverse direction of the tree", new ArrayList<>());
    
	
	
	final public Input<Double> scaleFactorInput = new Input<>("scaleFactor", "scaling factor: range from 0 to 1. Close to zero is very large jumps, close to 1.0 is very small jumps.", 0.75);
    final public Input<Boolean> rootOnlyInput = new Input<>("rootOnly", "scale a single node only", false);
    final public Input<Boolean> optimiseInput = new Input<>("optimise", "flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", true);
    final public Input<Double> scaleUpperLimit = new Input<>("upper", "Upper Limit of scale factor", 10.0);
    final public Input<Double> scaleLowerLimit = new Input<>("lower", "Lower limit of scale factor", 1e-8);
    final public Input<KernelDistribution> kernelDistributionInput = new Input<>("kernelDistribution", "provides sample distribution for proposals", KernelDistribution.newDefaultKernelDistribution());

    private double scaleFactor;
    protected KernelDistribution kernelDistribution;
    private double upper, lower;

    @Override
    public void initAndValidate() {
    	kernelDistribution = kernelDistributionInput.get();
        scaleFactor = scaleFactorInput.get();
        upper = scaleUpperLimit.get();
        lower = scaleLowerLimit.get();
    }


    protected double getScaler(int i, double value) {
		return kernelDistribution.getScaler(i, value, getCoercableParameterValue());
	}

    /**
     * override this for proposals,
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {

        try {

            double hastingsRatio;
            Stubs stubs = stubsInput.get();
            
            // Cache branch lengths before making proposal
            double[] cachedBranchLengths = stubs.prepareJacobian();
            
            final Tree tree = (Tree)InputUtil.get(treeInput, this);
            if (rootOnlyInput.get()) {
            	
            	
            	// Scale one node only
            	int nnodes = tree.getInternalNodeCount();
            	int nodeNr = Randomizer.nextInt(nnodes) + nnodes+1;
                //final Node node = tree.getNode(nodeNr);      
                final Node root = tree.getRoot();   
                if (root.isFake()) {
                	return Double.NEGATIVE_INFINITY; // Sampled ancestors
                }
                final double scale = getScaler(root.getNr(), root.getHeight());
                final double newHeight = root.getHeight() * scale;
                double parentHeight = root.isRoot() ? Double.POSITIVE_INFINITY : root.getParent().getHeight();
                if ( newHeight >= parentHeight || newHeight < Math.max(root.getLeft().getHeight(), root.getRight().getHeight())) {
                    return Double.NEGATIVE_INFINITY;
                }
                
                root.setHeight(newHeight);
                
                // Hastings ratio needs to be aware of how many are on the two branches under and one above 
                //int nstubs = stubs.getNStubsOnBranch(node.getLeft().getNr()) + stubs.getNStubsOnBranch(node.getRight().getNr());
                double stubJacobian = stubs.getLogJacobian(cachedBranchLengths);
               
                
                hastingsRatio = Math.log(scale) + stubJacobian;
                
            } else {
            	
            	final double scale = getScaler(0, Double.NaN);
            	
            	
            	// Keep track of branch length changes for hastings ratio
            	double[] oldLengths = new double[tree.getNodeCount()-1];
            	for (int i = 0; i < oldLengths.length; i ++) {
            		oldLengths[i] = tree.getNode(i).getLength();
            	}
            	
                // Scale the whole tree and a bunch of other parameters
                final int scaledNodes = tree.scale(scale);
                
                // Stubs are automatically scaled with the tree, as their heights are relative but hastings ratio needs to be aware
                double stubJacobian = stubs.getLogJacobian(cachedBranchLengths);
               
                
                int goingUp = scaledNodes;
                int goingDown = 0;
                for (StateNode up : upInput.get()) {
                    up = up.getCurrentEditable(this);
                    goingUp += up.scale(scale);
                }
               
                for (StateNode up : upInput.get()) {
                    if (outsideBounds(up)) {
                        return Double.NEGATIVE_INFINITY;
                    }
                }

                for (StateNode down : downInput.get()) {
                    down = down.getCurrentEditable(this);
                    goingDown += down.scale(1.0 / scale);
                }
                for (StateNode down : downInput.get()) {
                    if (outsideBounds(down)) {
                        return Double.NEGATIVE_INFINITY;
                    }
                }
                
                
                // Hastings ratio
                hastingsRatio = Math.log(scale) * (goingUp - goingDown) + stubJacobian;
               
                
            }


            return hastingsRatio;

        } catch (Exception e) {
            // whatever went wrong, we want to abort this operation...
            return Double.NEGATIVE_INFINITY;
        }
    }


    /**
     * automatic parameter tuning *
     */
    @Override
    public void optimize(final double logAlpha) {
        if (optimiseInput.get()) {
            double delta = calcDelta(logAlpha);
            delta += Math.log(1.0 / scaleFactor - 1.0);
            setCoercableParameterValue(1.0 / (Math.exp(delta) + 1.0));
        }
    }

    @Override
    public double getCoercableParameterValue() {
        return scaleFactor;
    }

    @Override
    public void setCoercableParameterValue(final double value) {
        scaleFactor = Math.max(Math.min(value, upper), lower);
    }

    @Override
    public String getPerformanceSuggestion() {
        final double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        final double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        final double sf = Math.pow(scaleFactor, ratio);

        final DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else if (prob > 0.40) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else return "";
    }
    
    
    private boolean outsideBounds(final StateNode node) {
        if (node instanceof Parameter<?>) {
            final Parameter<?> p = (Parameter<?>) node;
            final Double lower = (Double) p.getLower();
            final Double upper = (Double) p.getUpper();
            final Double value = (Double) p.getValue();
            if (value < lower || value > upper) {
                return true;
            }
        }
        return false;
    }
    
	
	@Override
    public List<StateNode> listStateNodes() {
		Stubs stubs = stubsInput.get();
        final List<StateNode> list = super.listStateNodes();
        list.add(stubs.getStubHeights());
        //list.add(stubs.getBranches());
        return list;
    }
	

} 

