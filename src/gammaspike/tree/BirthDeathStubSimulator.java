package gammaspike.tree;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.speciation.YuleModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.State;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;



@Description("Simulates a birth-death tree with psi sampling and labels the nodes. Assuming that rho=1")
public class BirthDeathStubSimulator extends YuleModel {
	
	final public Input<Integer> ntaxaInput = new Input<>("ntaxa", "number of leaves", Input.Validate.REQUIRED);
	final public Input<RealParameter> r0Input = new Input<>("r0", "ratio between birth and death rate", Input.Validate.OPTIONAL);
	final public Input<RealParameter> samplingProportionInput = new Input<>("samplingProportion", "sampling rate psi divided by (psi + mu)", Input.Validate.REQUIRED);
	final public Input<Boolean> removeStubsInput = new Input<>("removeStubs", "remove stubs from tree after simulation?", true);
	final public Input<Boolean> countOriginStubsInput = new Input<>("countOriginStubs", "include stubs on origin when counting stubs?", false);
	final public Input<IntegerParameter> stubsPerBranchInput = new Input<>("nstubsPerBranch", "number of stubs per branch (for simulation)", Input.Validate.OPTIONAL);
	
	
	final public static String NSTUBS_STR = "nstubs.branch";
	
	Tree fullTree;
	Tree reducedTree;
	int[] numberOfStubsPerBranch;
	int numberOfAncestralSamples = 0;
	
	
	@Override
    public void initAndValidate() {
		super.initAndValidate();
        fullTree = new Tree();
        reducedTree = new Tree();
	}
	
	
	/**
	 * Sample the tree
	 */
	private void sampleTree(State state, Random random) {
		
		
		// Standard Yule process
		if (r0Input.get() == null) {
        	super.sample(state, random);
        	return;
        }
        
		 
        // Get parameters and tree
        Tree tree = (Tree) treeInput.get();
        double birthRate = birthDiffRateParameterInput.get().getValue();
        double deathRate = birthRate / r0Input.get().getValue();
        double samplingProportion = samplingProportionInput.get().getValue();
        double samplingRate = samplingProportion*deathRate / (1 - samplingProportion);
        double pbirth = birthRate / (birthRate + deathRate);
        
        

        
        while (true) {
        	
        	
            
            // Sample until we have ntaxa
            int targetN = ntaxaInput.get();
            Node root = new Node();
            root.setHeight(0);
            double time = 0;
            List<Node> livingLeaves = new ArrayList<>();
            List<Node> extinctLeaves = new ArrayList<>();
            List<Node> allNodes = new ArrayList<>();
            livingLeaves.add(root);
            allNodes.add(root);
	        
	        while (livingLeaves.size() < targetN) {
	        	
	        	
	        	
	        	
	        	// Extinct - try again
	        	if (livingLeaves.isEmpty()) {
	        		//Log.warning("Gone extinct. Trying again.");
	        		break;
	        	}
	        	
	        	// Total rate
	        	double rate = livingLeaves.size() * (birthRate + deathRate);
	        	
	        	
	        	// Increment time
	        	double dt = Randomizer.nextExponential(rate);
	        	
	        	
	        	// Make sure all living leaves have the new time
	        	for (Node leaf : livingLeaves) {
	        		leaf.setHeight(dt + time);
	        	}
	        	
	        	// Birth or death?
	        	boolean birth = Randomizer.nextFloat() < pbirth;
	        	if (birth) {
	
	        		// Add 2 new nodes
	        		Node child1 = new Node();
	        		Node child2 = new Node();
	        		child1.setHeight(time + dt);
	        		child2.setHeight(time + dt);
	        		
	        		
	        		// Sample their parent
	        		int parentNr = Randomizer.nextInt(livingLeaves.size());
	        		Node parent = livingLeaves.get(parentNr);
	        		parent.addChild(child1);
	        		parent.addChild(child2);
	        		
	        		// Update leaf list
	        		livingLeaves.remove(parent);
	        		livingLeaves.add(child1);
	        		livingLeaves.add(child2);
	        		
	        		allNodes.add(child1);
	        		allNodes.add(child2);
	        		
	        		
	        		
	        	}else {
	        		
	        		// Remove a node from leaf list
	        		int deceasedNr = Randomizer.nextInt(livingLeaves.size());
	        		Node deceased = livingLeaves.get(deceasedNr);
	        		livingLeaves.remove(deceased);
	        		extinctLeaves.add(deceased);
	        		
	        	}
	        	
	        	
	        	
	        	time = time + dt;
	        	
	        	//Log.warning("time " + time + " dt " + dt);
	        	
	        	
	        }
	        
	        
	        // Try again
	        if (livingLeaves.isEmpty()) {
	        	continue;
	        }
	        

	        
	        // Add a little more onto each leaf
	        double rate = livingLeaves.size() * (birthRate + deathRate);
        	double dt = Randomizer.nextExponential(rate);
        	time += dt;
	        for (Node node : livingLeaves) {
	        	if (node.isLeaf()) node.setHeight(time);
	        }
	        
	        
	
	        // Reverse time
	        for (Node node : allNodes) {
	        	node.setHeight(time - node.getHeight());
	        }
	        
	        
	
	        // Number extant leaves 1-N
	        int count = 0;
	        for (Node node : livingLeaves) {
	    		node.setNr(count);
	    		//node.setID("taxon" + (count+1));
	    		count ++;
	        }
	        
	        // Number extinct leaves next
	        for (Node node : extinctLeaves) {
	    		node.setNr(count);
	    		//node.setID("taxon" + (count+1));
	    		count ++;
	        }
	        
	        
	        // Internal nodes next
	        for (Node node : allNodes) {
	        	if (!node.isLeaf() && !node.isRoot()) {
	        		node.setNr(count);
	        		count ++;
	        	}
	        }
	        
	        
	        // Root last
	        root.setNr(count);
	        
	        
	        
	        // If one of the two root children are extinct, try again
	        if (allDescendantsAreUnsampled(root.getLeft()) || allDescendantsAreUnsampled(root.getRight())) {
	        	Log.warning("bad root");
	        	continue;
	        }
	        
	        
	        
	        // Add samples onto old part of tree
	        numberOfAncestralSamples = 0;
	        this.sampleAncestralLineages(root, samplingRate);
	        
	        
	        

	        // Remove unsampled lineages
        	Tree prunedTree = this.removeStubs(root.copy());
        	
        	
        	 // Temp: set all sampled ancestors to standard extinct leaves
//	        for (Node node : prunedTree.getNodesAsArray()) {
//	        	if (node.isDirectAncestor()) {
//	        		double h = node.getHeight();
//	        		node.setHeight(h * 0.95);
//	        	}
//	        }
	        
        	
        	
        	//Log.warning("total nleaves after pruning unsampled lineages: " + prunedTree.getRoot().getLeafNodeCount());
        	reducedTree.assignFrom(prunedTree);
        	Tree newTree = new Tree(root);
	        fullTree.assignFrom(newTree);
	        
	        
	        
	       
	        
	        
	        
	        if (removeStubsInput.get()) {
	        	
	        	// Remove extinct lineages
	        	tree.assignFrom(prunedTree);     
	        	
	        }else {
	        	
	        	// Keep extinct lineages
		        tree.assignFrom(newTree);
	        }
	        
	        
	        // New number of leaves
	        int nleaves = tree.getLeafNodeCount();

	        
	        
	        // Assign names
	        //for (Node node : prunedTree.getRoot().getAllChildNodesAndSelf()) {
	        	//if (node.isLeaf()) {
	        		//node.setID("ataxon" + node.getNr());
	        	//}
	       // }
	        

	        
	        // Stub counting
	        if (stubsPerBranchInput.get() != null) {
	        	stubsPerBranchInput.get().setDimension(prunedTree.getNodeCount());
	        }
	        numberOfStubsPerBranch = new int[prunedTree.getNodeCount()];
	        for (int i = 0; i < numberOfStubsPerBranch.length;  i++) {
	        	
	        	int nstubs = (int)(prunedTree.getNode(i).getMetaData(NSTUBS_STR));
	        	numberOfStubsPerBranch[i] = nstubs; 
	        	if (stubsPerBranchInput.get() != null) {
	        		stubsPerBranchInput.get().setValue(i, nstubs);
	        	}
	        	
	        	
	        }
	        
	        
	        //Log.warning("Success." + tree);
	        return;
        }
        
	}
	
	
	/**
	 * Add sampled ancestors onto branches (i.e., leaves with branch length 0)
	 * @param node
	 * @param samplingRate
	 */
	private void sampleAncestralLineages(Node node, double samplingRate) {
		
		
		if (samplingRate <= 0) return;
		
		// Take the children now before they change
		Node leftChild = null, rightChild = null;
		if (!node.isLeaf()) {
			leftChild = node.getChild(0);
			rightChild = node.getChild(1);
		}
		
		
		if (!node.isRoot()) {
			
			// How many samples on this branch?
			double length = node.getLength();
			double height = node.getHeight();
			double poissonRate = length*samplingRate;
			int numberOfSamples = (int)Randomizer.nextPoisson(poissonRate);
			numberOfAncestralSamples += numberOfSamples;
			//Log.warning("Adding " + numberOfSamples + " samples ");
			if (numberOfSamples > 0) {
				
				
				// Sample times
				List<Double> timesOfSamples = new ArrayList<>();
				for (int i = 0; i < numberOfSamples; i++) {
					double time = Randomizer.nextDouble()*length + height;
					timesOfSamples.add(time);
				}
				Collections.sort(timesOfSamples);
				
				// Create new node - going up the branch from young to old
				for (int i = 0; i < numberOfSamples; i++) {
					
					
					double time = timesOfSamples.get(i);
					
					if (time < node.getHeight()) {
						throw new IllegalArgumentException("Dev error 3535: negative branch length detected");
					}
					
					// Make leaf
					Node sampledLeaf = new Node();
					sampledLeaf.setHeight(time);
					

					
					// Make new internal node
					Node internal = new Node();
					internal.setHeight(time);
					internal.addChild(sampledLeaf);
					
					
					// Add to the branch
					Node parent = node.getParent();
					//Log.warning("nleaves before " + parent.getLeafNodeCount() + " time " + time);
					
					
					// Set parent
					parent.removeChild(node);
					node.setParent(null);
					internal.addChild(node);
					parent.addChild(internal);
					node = internal;
					//Log.warning("nleaves after " + parent.getLeafNodeCount());
					
					if (node.getLength() < 0) {
						throw new IllegalArgumentException("Dev error 3536: negative branch length detected");
					}
					
					
				}
				
				
				
			}
			
			
		}
		
		
		if (leftChild != null) {
			sampleAncestralLineages(leftChild, samplingRate);
			sampleAncestralLineages(rightChild, samplingRate);
		}
		
		
		if (node.isRoot()) {
			
			// Renumber
			int nodeNr = 0;
			for (Node leaf : node.getAllLeafNodes()) {
				leaf.setNr(nodeNr);
				nodeNr++;
			}
			for (Node internal : node.getAllChildNodesAndSelf()) {
				if (internal.isRoot() || internal.isLeaf()) continue;
				internal.setNr(nodeNr);
				nodeNr++;
			}
			node.setNr(nodeNr);
			
		}
		
		
		
		
	}
	
	
	// Create a new tree with stubs removed
	private Tree removeStubs(Node root) {
		
		//System.out.println("nc1=" + root.getNodeCount());
		// Set the root such that both children are extant
		while (true) {
			
			if (root.isLeaf()) break;
		
			Node left = root.getChild(0);
			Node right = root.getChild(1);
			if (allDescendantsAreUnsampled(left) && allDescendantsAreUnsampled(right)) {
				throw new IllegalArgumentException("Unexpected error: all children are extinct");
			}
			
			if (allDescendantsAreUnsampled(left)) {
				root = right;
			}else if (allDescendantsAreUnsampled(right)) {
				root = left;
			}else {
				break;
			}
		
		}
		root.setParent(null);
		//System.out.println("nc2=" + root.getNodeCount());
		removeStubsRecursive(root, 0);
		
		// Renumber
		int nodeNr = 0;
		for (Node leaf : root.getAllLeafNodes()) {
			leaf.setNr(nodeNr);
			nodeNr++;
		}
		for (Node internal : root.getAllChildNodesAndSelf()) {
			if (internal.isRoot() || internal.isLeaf()) continue;
			internal.setNr(nodeNr);
			nodeNr++;
		}
		root.setNr(nodeNr);
		
		Tree reducedTree = new Tree(root);
		return reducedTree;
	}
	
	public void removeStubsRecursive(Node node, int nStubsOnBranch) {
		
		if (node.isLeaf()) {
			node.setMetaData(NSTUBS_STR, nStubsOnBranch);
			return;
		}
		
		Node parent = node.getParent();
		Node left = node.getChild(0);
		Node right = node.getChild(1);
		
		
		// Should never happen by the design of this recursion
		if (allDescendantsAreUnsampled(left) && allDescendantsAreUnsampled(right)) {
			throw new IllegalArgumentException("Unexpected error: all children are extinct");
		}
		
		
		// Append the right child to this parent, and remove this node
		if (allDescendantsAreUnsampled(left)) {
			parent.removeChild(node);
			node.setParent(null);
			parent.addChild(right);
			removeStubsRecursive(right, nStubsOnBranch+1);
		}
		
		// Append the left child to this parent, and remove this node
		else if (allDescendantsAreUnsampled(right)) {
			parent.removeChild(node);
			node.setParent(null);
			parent.addChild(left);
			removeStubsRecursive(left, nStubsOnBranch+1);
		}
		
		
		// Recurse along both children
		else {
			
			node.setMetaData(NSTUBS_STR, nStubsOnBranch);
			removeStubsRecursive(left, 0);
			removeStubsRecursive(right, 0);
		}
		
	}


	@Override
    public void sample(State state, Random random) {

        if (sampledFlag) return;
        
        
        // Sample the tree
        sampleConditions(state, random);
        sampleTree(state, random);
        sampledFlag = true;
       
        
	}
	
	
	
	@Override
    public void init(final PrintStream out) {
		out.print("nancestors\t");
        out.print("nstubs\t");
        for (int i = 0; i < ntaxaInput.get()*2-1; i ++) {
        	out.print("nstubs.branch." + (i+1) + "\t");
        }
    }

    @Override
    public void log(final long sample, final PrintStream out) {
    	
    	// Number of ancestral samples
    	out.print(numberOfAncestralSamples + "\t");
        
        // Number of stubs
        out.print(getNstubs() + "\t");

        
        // Number of bursts per branch
        for (int i = 0; i < ntaxaInput.get()*2-1; i ++) {
        	out.print(numberOfStubsPerBranch[i] + "\t");
        }
        
    }
    
    
    public int getNstubs() {
    	
    	Tree tree = fullTree;
    	Node root = tree.getRoot();
    	if (!countOriginStubsInput.get()) {
    		root = getEmpiricalRoot(root);
    	}
    	int nstubs = this.countStubs(root, 0);
		return nstubs;
	}

    
    
    /**
     * Find the oldest node in the tree such that both child lineages are non-extinct
     * @param node
     * @return
     */
    private Node getEmpiricalRoot(Node node) {
    	
    	
    	if (node.isLeaf()) return node;
    	
    	Node left = node.getChild(0);
		Node right = node.getChild(1);
		
		
		// Found the root
		if (!allDescendantsAreUnsampled(left) && !allDescendantsAreUnsampled(right)) {
			return node;
		}
		
		
		// Left is extinct, so root is on right side
		if (allDescendantsAreUnsampled(left)) {
			return getEmpiricalRoot(right);
		}
		
		// Right is extinct, so root is on left side
		else if (allDescendantsAreUnsampled(right)) {
			return getEmpiricalRoot(left);
		}
		
		
		throw new IllegalArgumentException("Unexpected error: all children are extinct");
    	
    	
    }

	private int countStubs(Node node, int nstubs) {
		
		
		if (allDescendantsAreUnsampled(node)) {
			return nstubs + 1;
		}
		
		if (node.isLeaf()) {
			return nstubs;
		}
		
		
		for (Node child : node.getChildren()) {
			nstubs = countStubs(child, nstubs);
		}
		
		
		return nstubs;
	}


	
	/**
	 * A node is unsampled if a) it is not extant (height>0), and b) it has branch length > 0 
	 * @param node
	 * @return
	 */
	private boolean allDescendantsAreUnsampled(Node node) {
		
		
		if (node.isLeaf()) {
			if (node.getHeight() > 0 && node.getLength() > 1e-16) return true;
			return false;
		}
		
		
		boolean allChildrenExtinct = true;
		for (Node child : node.getChildren()) {
			if (!allDescendantsAreUnsampled(child)){
				allChildrenExtinct = false;
				break;
			}
		}
		
		return allChildrenExtinct;
		
	}


	@Override
    public List<String> getConditions() {
        List<String> conditions = new ArrayList<>();
        conditions.add(birthDiffRateParameterInput.get().getID());
        conditions.add(samplingProportionInput.get().getID());
        conditions.add(r0Input.get().getID());
        return conditions;
    }


	

}



