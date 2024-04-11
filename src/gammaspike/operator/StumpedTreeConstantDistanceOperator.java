package gammaspike.operator;


import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.operator.kernel.KernelDistribution;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import gammaspike.tree.Stubs;


@Description("InConstantDistanceOperator from the ORC package but with stub awareness")
public class StumpedTreeConstantDistanceOperator extends TreeOperator {

	final public Input<Stubs> stubsInput = new Input<>("stubs", "stubs for the tree", Input.Validate.REQUIRED);
	final public Input<Double> twindowSizeInput = new Input<>("twindowSize", "the size of the window when proposing new node time", Input.Validate.REQUIRED);
	final public Input<RealParameter> rateInput = new Input<>("rates", "the rates associated with nodes in the tree for sampling of individual rates among branches.", Input.Validate.REQUIRED);
	final public Input<KernelDistribution> proposalKernelInput = new Input<>("kernel", "Proposal kernel for a random walk on the internal node height.", KernelDistribution.newDefaultKernelDistribution());
	
	// Proposal kernel
	private KernelDistribution kernel;
	
	private double twindowSize;
	private RealParameter rates;

	
	@Override
    public void initAndValidate() {
    	this.twindowSize = twindowSizeInput.get();
    	this.rates = rateInput.get();
        this.kernel = proposalKernelInput.get();
    }
	
	
	@Override
    public double proposal() {
		
		// Cache branch lengths before making proposal
		Stubs stubs = stubsInput.get();
        double[] cachedBranchLengths = stubs.prepareJacobian();
		
		double HR = proposalInner();
		
		// Jacobian. Relative stub heights stay the same but absolute heights change
        double logJacobian = stubs.getLogJacobian(cachedBranchLengths);
        
        return HR + logJacobian;
		
	}
	
	
    public double proposalInner() {
    	
        final Tree tree = treeInput.get();
        int nodeCount = tree.getNodeCount(); //return the number of nodes in the tree
        int branchCount = nodeCount - 1; //the number of branches of the tree

        // the chosen node to work on
        Node node;

        // The original node times
        double t_x, t_j, t_k;

        // The original rates
        double r_x, r_k, r_j;


        double hastingsRatio = 0.0;

        // Randomly select an internal node, denoted by node x avoiding fake nodes used to connect direct ancestors into tree.
        node = this.sampleNode(tree);
        if (node == null) return Double.NEGATIVE_INFINITY;

       // the number of this node
        int nodeNr = node.getNr();
        // if this node has max number, then use the free index stored in root node to get rate.
        if (nodeNr == branchCount) {
            nodeNr = node.getTree().getRoot().getNr();
        }

       // time for this node
       t_x = node.getHeight();


       // Step 2: Access to the child nodes of this node
       // son
       Node son = node.getChild(0); // get the left child of this node, i.e. son
       t_j = son.getHeight(); // node time of son
       int sonNr = son.getNr(); // node number of son
       if (sonNr == branchCount) {
           sonNr = son.getTree().getRoot().getNr();
        }

       // daughter
       Node daughter = node.getChild(1); //get the right child of this node, i.e. daughter
       t_k = daughter.getHeight(); // node time of daughter
       int dauNr = daughter.getNr(); // node number of daughter
       if (dauNr == branchCount) {
            dauNr = daughter.getTree().getRoot().getNr();
       }

       r_x = rates.getValue(nodeNr); // rate of branch above this node
       r_j = rates.getValue(sonNr); // rate of branch above son
       r_k = rates.getValue(dauNr); // rate of branch above daughter


       // Propose a new node time 
       double a;
       if (kernel != null) a = kernel.getRandomDelta(1, twindowSize);
       else a = Randomizer.uniform(-twindowSize, twindowSize);
       double t_x_ = t_x + a;

       // deal with the boundary cases
       double upper = node.getParent().getHeight();
       double lower = Math.max(t_j, t_k);

        if (t_x_ == lower || t_x_ == upper) {
            return Double.NEGATIVE_INFINITY;
        }
        // fold the proposed node time
        double err; double n; double r;
        if (t_x_ > upper) {
            err = t_x_ - upper;
            n = Math.floor(err / (upper - lower));
            r = err - n * (upper - lower);
            if (n % 2 == 0) {
                t_x_ = upper - r;
            } else {
                t_x_ = lower + r;
            }
        } else if (t_x_ < lower) {
            err = lower - t_x_;
            n = Math.floor(err / (upper - lower));
            r = err - n * (upper - lower);
            if (n % 2 == 0) {
                t_x_ = lower + r;
            } else {
                t_x_ = upper - r;
            }
        }

        // set the proposed node time
        node.setHeight(t_x_);


       // Propose the new rates
       // r_x_, r_j_, r_k_
       double r_x_ = r_x * (upper - t_x) / (upper - t_x_);
       double r_j_ = r_j * (t_x - t_j) / (t_x_ - t_j);
       double r_k_ = r_k * (t_x - t_k) / (t_x_ - t_k);


       // Set the proposed new rates 
       rates.setValue(nodeNr, r_x_);
       rates.setValue(sonNr, r_j_);
       rates.setValue(dauNr, r_k_);

        

        // Calculate Hastings ratio
        double nu =(upper - t_x) * (t_x - t_j) * (t_x - t_k) ;
        double de = (upper - t_x_) * (t_x_ - t_j) * (t_x_ - t_k);
        hastingsRatio = Math.log(nu / de);
        return hastingsRatio;
        
    }
    
    

    protected Node sampleNode(Tree tree) {
    	Node node = null;
    	int nodeCount = tree.getNodeCount();
    	
    	if (tree.getLeafNodeCount() < 3) return null;
    	
    	do {
            final int nodeNr = nodeCount / 2 + 1 + Randomizer.nextInt(nodeCount / 2);
            node = tree.getNode(nodeNr);
    	} while (node.isRoot() || node.isLeaf() || node.isFake());
    	
    	
    	return node;
    	 
	}


	
}





