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
		
		// If no r0 is provided, use the standard Yule process
		if (r0Input.get() == null) {
        	super.sample(state, random);
        	return;
        }

        // Get parameters and tree
        Tree tree = (Tree) treeInput.get();
        double birthRate = birthDiffRateParameterInput.get().getValue();
        double deathRate = birthRate / r0Input.get().getValue();
        double samplingProportion = samplingProportionInput.get().getValue();
        double samplingRate = samplingProportion * deathRate / (1 - samplingProportion);
        double probBirth = birthRate / (birthRate + deathRate);

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
	        	boolean birth = Randomizer.nextFloat() < probBirth;
	        	if (birth) { // Node splits into two children

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

	        	} else { // Node dies

	        		// Remove a node from leaf list
	        		int deceasedNr = Randomizer.nextInt(livingLeaves.size());
	        		Node deceased = livingLeaves.get(deceasedNr);
	        		livingLeaves.remove(deceased);
	        		extinctLeaves.add(deceased);
	        		
	        	}
	        	
	        	time = time + dt;
	        	//Log.warning("time " + time + " dt " + dt);
	        	
	        }


	        // Restart the simulation if all lineages have gone extinct
	        if (livingLeaves.isEmpty()) {
	        	continue;
	        }


	        // Add a little more time onto each leaf node
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
	        

	        // Number extant leaves
	        int count = 0;
	        for (Node node : livingLeaves) {
	    		node.setNr(count);
	    		//node.setID("taxon" + (count+1));
	    		count ++;
	        }
	        // Number extinct leaves
	        for (Node node : extinctLeaves) {
	    		node.setNr(count);
	    		//node.setID("taxon" + (count+1));
	    		count ++;
	        }
	        // Number internal nodes
	        for (Node node : allNodes) {
	        	if (!node.isLeaf() && !node.isRoot()) {
	        		node.setNr(count);
	        		count ++;
	        	}
	        }
	        // Number root node
	        root.setNr(count);
	        
	        
			// Restart the simulation if one of the two root children has no sampled descendants
	        if (allDescendantsAreUnsampled(root.getLeft()) || allDescendantsAreUnsampled(root.getRight())) {
	        	//Log.warning("bad root");
	        	continue;
	        }
	        
	        
			// Simulate sampling through time
	        numberOfAncestralSamples = 0;
	        this.sampleAncestralLineages(root, samplingRate);
	        
	        
	        // Remove unsampled lineages
        	Tree prunedTree = this.removeStubs(root.copy());
        	

        	reducedTree.assignFrom(prunedTree); // reducedTree: Pruned tree with only sampled lineages
        	Tree newTree = new Tree(root);
	        fullTree.assignFrom(newTree); // fullTree: Tree with all sampled and unsampled nodes
	        
	        
	        if (removeStubsInput.get()) {
	        	// Remove unsampled lineages
	        	tree.assignFrom(prunedTree);     
	        } else {
	        	// Keep unsampled lineages
		        tree.assignFrom(newTree);
	        }
	        

	        // Stub counting
	        if (stubsPerBranchInput.get() != null) {
	        	stubsPerBranchInput.get().setDimension(prunedTree.getNodeCount());
	        }
	        numberOfStubsPerBranch = new int[prunedTree.getNodeCount()];
			// Loop over each node in the pruned tree
	        for (int i = 0; i < numberOfStubsPerBranch.length;  i++) {
				// Retrieves the stub count stored as metadata (NSTUBS_STR) and stores it in the array
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
	 * If sampled ancestors have no sampled descendants, they become sampled tips through time
	 * @param node
	 * @param samplingRate
	 */
	private void sampleAncestralLineages(Node node, double samplingRate) {

		// If the sampling rate is 0 or negative, do nothing (i.e., no ancestral samples).
		if (samplingRate <= 0) return;

		// Save the children of the current node before the tree structure is modified
		Node leftChild = null, rightChild = null;
		if (!node.isLeaf()) {
			leftChild = node.getChild(0);
			rightChild = node.getChild(1);
		}

		// Sample along branch
		if (!node.isRoot()) {
			
			// Calculate number of ancestral samples
			double length = node.getLength(); // Branch length
			double height = node.getHeight(); // Node height
			double poissonRate = length * samplingRate;
			int numberOfSamples = (int)Randomizer.nextPoisson(poissonRate); // Ancestor sampling is Poisson-distributed along the branch
			
			numberOfAncestralSamples += numberOfSamples;

			if (numberOfSamples > 0) {

				// Sample each sampled ancestor at a random time along the branch
				List<Double> timesOfSamples = new ArrayList<>();
				for (int i = 0; i < numberOfSamples; i++) {
					double time = Randomizer.nextDouble() * length + height;
					timesOfSamples.add(time);
				}
				Collections.sort(timesOfSamples);
				
				// Insert a new sampled ancestor at each sampled time
				for (int i = 0; i < numberOfSamples; i++) {
					double time = timesOfSamples.get(i);
					
					if (time < node.getHeight()) {
						throw new IllegalArgumentException("Dev error 3535: negative branch length detected");
					}
					
					// Create a new leaf node at the sampled time
					Node sampledLeaf = new Node();
					sampledLeaf.setHeight(time);
					
					// Create a new internal node at the sampled time, act as a branching point between the sampled leaf and the tree
					Node internal = new Node();
					internal.setHeight(time);
					internal.addChild(sampledLeaf);

					// Insert the sampled ancestor to the tree
					Node parent = node.getParent(); // Get the parent node of the current node
					parent.removeChild(node); // Temporarily detach the current node
					node.setParent(null);
					internal.addChild(node); // Make the current node a child of the internal node
					parent.addChild(internal); // Make the internal node a child of the original parent node: parent → internal → node
					node = internal; // Reassign node to the new internal node

					if (node.getLength() < 0) {
						throw new IllegalArgumentException("Dev error 3536: negative branch length detected");
					}
					
				}
				
			}

		}

		// Recursively applies the same sampling process down both subtrees.
		if (leftChild != null) {
			sampleAncestralLineages(leftChild, samplingRate);
			sampleAncestralLineages(rightChild, samplingRate);
		}

		// Renumber nodes
		if (node.isRoot()) {
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
		
		// Set the root such that both child nodes are extant/have sampled descendants
		while (true) {
			
			if (root.isLeaf()) break;
		
			Node left = root.getChild(0);
			Node right = root.getChild(1);
			// Throws an error if both child branches are unsampled
			if (allDescendantsAreUnsampled(left) && allDescendantsAreUnsampled(right)) {
				throw new IllegalArgumentException("Unexpected error: all children are extinct");
			}
			if (allDescendantsAreUnsampled(left)) {
				root = right;
			} else if (allDescendantsAreUnsampled(right)) {
				root = left;
			} else {
				break;
			}

		}

		root.setParent(null);
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

		// Throws an error if both child branches have no sampled descendant
		if (allDescendantsAreUnsampled(left) && allDescendantsAreUnsampled(right)) {
			throw new IllegalArgumentException("Unexpected error: all children are extinct");
		}
		// If left child has no sampled descendant
		if (allDescendantsAreUnsampled(left)) {
			parent.removeChild(node);
			node.setParent(null);
			parent.addChild(right); // Replace with the sampled right child
			removeStubsRecursive(right, nStubsOnBranch+1);
		}
		// If right child has no sampled descendant
		else if (allDescendantsAreUnsampled(right)) {
			parent.removeChild(node);
			node.setParent(null);
			parent.addChild(left); // Replace with the sampled left child
			removeStubsRecursive(left, nStubsOnBranch+1);
		}
		// If neither child is unsampled, store the number of stubs, recurse down both children
		else {
			node.setMetaData(NSTUBS_STR, nStubsOnBranch);
			removeStubsRecursive(left, 0);
			removeStubsRecursive(right, 0);
		}

	}


	/**
	 * A node is unsampled if
	 * a) it is not extant (height > 0), and
	 * b) it has branch length > 0
	 * or, c) all children are extinct without samples
	 * @param node
	 * @return
	 */
	private boolean allDescendantsAreUnsampled(Node node) {

		if (node.isLeaf()) {
			if (node.getHeight() > 0 && node.getLength() > 1e-16) return true;
			return false;
		}

		boolean allChildrenExtinct = true;
		// Loop over all children
		for (Node child : node.getChildren()) {
			// If any child has sampled descendants
			if (!allDescendantsAreUnsampled(child)){
				allChildrenExtinct = false;
				break;
			}
		}
		return allChildrenExtinct;

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



	@Override
    public List<String> getConditions() {
        List<String> conditions = new ArrayList<>();
        conditions.add(birthDiffRateParameterInput.get().getID());
        conditions.add(samplingProportionInput.get().getID());
        conditions.add(r0Input.get().getID());
        return conditions;
    }


	

}
