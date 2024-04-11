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



@Description("Simulates a birth-death tree and labels the nodes")
public class ForwardTimeSimulatorResub extends YuleModel {
	
	final public Input<Integer> ntaxaInput = new Input<>("ntaxa", "number of leaves", Input.Validate.REQUIRED);
	//final public Input<RealParameter> deathRateParameterInput = new Input<>("mu", "death rate", Input.Validate.OPTIONAL);
	final public Input<RealParameter> r0Input = new Input<>("r0", "ratio between birth and death rate", Input.Validate.OPTIONAL);
	final public Input<RealParameter> pInput = new Input<>("p", "probability one child is the same as parent", Input.Validate.REQUIRED);
	final public Input<IntegerParameter> nInput = new Input<>("N", "total space of alphabet", Input.Validate.REQUIRED);
	final public Input<IntegerParameter> statesInput = new Input<>("states", "the state of each node in the tree", Input.Validate.REQUIRED);
	
	final public Input<Boolean> removeStumpsInput = new Input<>("removeStumps", "remove stumps from tree after simulation?", true);
	final public Input<Boolean> countOriginStubsInput = new Input<>("countOriginStubs", "include stubs on origin when counting stubs?", false);
	final public Input<IntegerParameter> stubsPerBranchInput = new Input<>("nstubsPerBranch", "number of stubs per branch (for simulation)", Input.Validate.OPTIONAL);
	
	
	List<Integer> leafStates = new ArrayList<>();
	int[] leafDistancesFromRoot;
	boolean[] leafReps;
	
	final public static String NSTUBS_STR = "nstubs.branch";
	
	Tree fullTree;
	Tree reducedTree;
	int[] numberOfStubsPerBranch;
	
	
	
	@Override
    public void initAndValidate() {
		super.initAndValidate();
		leafStates = new ArrayList<>();
        leafDistancesFromRoot = new int[ntaxaInput.get()];
        leafReps = new boolean[ntaxaInput.get()];
        fullTree = new Tree();
        reducedTree = new Tree();
	}
	
	
	/**
	 * Sample the tree
	 */
	private void sampleTree(State state, Random random) {
		
		
		// Standatd Yule process
		if (r0Input.get() == null) {
        	super.sample(state, random);
        	return;
        }
        
		 
        // Get parameters and tree
        Tree tree = (Tree) treeInput.get();
        double birthRate = birthDiffRateParameterInput.get().getValue();
        double deathRate = birthRate / r0Input.get().getValue();
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
	    		node.setID("taxon" + (count+1));
	    		count ++;
	        }
	        
	        // Number extinct leaves next
	        for (Node node : extinctLeaves) {
	    		node.setNr(count);
	    		node.setID("taxon" + (count+1));
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
	        if (allDescendantsAreExtinct(root.getLeft()) || allDescendantsAreExtinct(root.getRight())) {
	        	Log.warning("bad root");
	        	continue;
	        }
	        
	        

	        // Remove extinct lineages
        	Tree prunedTree = this.removeStumps(root.copy());
        	reducedTree.assignFrom(prunedTree);
        	Tree newTree = new Tree(root);
	        fullTree.assignFrom(newTree);
	        
	        if (removeStumpsInput.get()) {
	        	
	        	// Remove extinct lineages
	        	tree.assignFrom(prunedTree);     
	        	
	        }else {
	        	
	        	// Keep extinct lineages
		        tree.assignFrom(newTree);
	        }
	        
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
	
	
	
	// Create a new tree with stumps removed
	private Tree removeStumps(Node root) {
		
		//System.out.println("nc1=" + root.getNodeCount());
		// Set the root such that both children are extant
		while (true) {
			
			if (root.isLeaf()) break;
		
			Node left = root.getChild(0);
			Node right = root.getChild(1);
			if (allDescendantsAreExtinct(left) && allDescendantsAreExtinct(right)) {
				throw new IllegalArgumentException("Unexpected error: all children are extinct");
			}
			
			if (allDescendantsAreExtinct(left)) {
				root = right;
			}else if (allDescendantsAreExtinct(right)) {
				root = left;
			}else {
				break;
			}
		
		}
		root.setParent(null);
		//System.out.println("nc2=" + root.getNodeCount());
		removeStumpsRecursive(root, 0);
		
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
	
	public void removeStumpsRecursive(Node node, int nStubsOnBranch) {
		
		if (node.isLeaf()) {
			node.setMetaData(NSTUBS_STR, nStubsOnBranch);
			return;
		}
		
		Node parent = node.getParent();
		Node left = node.getChild(0);
		Node right = node.getChild(1);
		
		
		// Should never happen by the design of this recursion
		if (allDescendantsAreExtinct(left) && allDescendantsAreExtinct(right)) {
			throw new IllegalArgumentException("Unexpected error: all children are extinct");
		}
		
		
		// Append the right child to this parent, and remove this node
		if (allDescendantsAreExtinct(left)) {
			parent.removeChild(node);
			node.setParent(null);
			parent.addChild(right);
			removeStumpsRecursive(right, nStubsOnBranch+1);
		}
		
		// Append the left child to this parent, and remove this node
		else if (allDescendantsAreExtinct(right)) {
			parent.removeChild(node);
			node.setParent(null);
			parent.addChild(left);
			removeStumpsRecursive(left, nStubsOnBranch+1);
		}
		
		
		// Recurse along both children
		else {
			
			node.setMetaData(NSTUBS_STR, nStubsOnBranch);
			removeStumpsRecursive(left, 0);
			removeStumpsRecursive(right, 0);
		}
		
	}


	@Override
    public void sample(State state, Random random) {

        if (sampledFlag) return;
        
        
        // Sample the tree
        sampleConditions(state, random);
        sampleTree(state, random);
        sampledFlag = true;
       
        
        Tree tree = fullTree;
        leafStates.clear();

        
        
        int n = nInput.get().getValue();
        double p = pInput.get().getValue();
        IntegerParameter states = statesInput.get();
        states.setDimension(tree.getNodeCount());
        
        
        // Simulate forward time labelling process from the root
        int stateNr = random.nextInt(n);
        Node node = tree.getRoot();
        states.setValue(node.getNr(), stateNr);
        
        
        int inheritState = random.nextFloat() < p ? stateNr : random.nextInt(n);
        sampleLabelsInClade(node.getLeft(), p, n, states, inheritState, random);
        sampleLabelsInClade(node.getRight(), p, n, states, -1, random);
        
        
        // Reverse time accumulation: what states are in the root?
        List<Integer> rootStates = getStatesInClade(node, states);
        

        
        // For each leaf, is its state represented in the root?
        for (int leafNr = 0; leafNr < tree.getLeafNodeCount(); leafNr++) {
        	if (tree.getNode(leafNr).getHeight() > 0) continue;
        	int leafState = states.getValue(leafNr);
        	boolean presentInRoot = rootStates.contains(leafState);
        	leafReps[leafNr] = presentInRoot;
        }
        
	}
	
	
	private List<Integer> getStatesInClade(Node node, IntegerParameter statesOnTree){
		
		List<Integer> relabelledStates = new ArrayList<>();
		if (node.isLeaf()) {
			 relabelledStates.add(statesOnTree.getValue(node.getNr()));
			 return relabelledStates;
		}
		
		
		List<Integer> leftStates = getStatesInClade(node.getLeft(), statesOnTree);
		List<Integer> rightStates = getStatesInClade(node.getRight(), statesOnTree);
		
		
		int stateNr = statesOnTree.getValue(node.getNr());
		int leftNr = statesOnTree.getValue(node.getLeft().getNr());
		int rightNr = statesOnTree.getValue(node.getRight().getNr());
		
		
		// Merge
		if (stateNr != leftNr && stateNr != rightNr) {
			relabelledStates.addAll(leftStates);
			for (int state : rightStates) {
				if (!relabelledStates.contains(state)) {
					relabelledStates.add(state);
				}
			}
			Collections.sort(relabelledStates);
			
			//Log.warning("r " + stateNr + " -> " + leftNr + " ," + rightNr);
			
		}else {
			
			
			// Just take left
			relabelledStates.addAll(leftStates);
			
			//Log.warning("exp " + stateNr + " -> " + leftNr + " ," + rightNr);
			
		}
		
		return relabelledStates;
		
	}
	
	
	private void sampleLabelsInClade(Node node, double p, int n, IntegerParameter states, int inheritState, Random random) {
		
		
		// Value of this state
		int stateNr = inheritState > -1 ? inheritState  : random.nextInt(n);
		states.setValue(node.getNr(), stateNr);
				
		
		if (node.isLeaf()) {
			if (!leafStates.contains(stateNr)) {
				leafStates.add(stateNr);
			}
			
			if (node.getHeight() > 0) return;
			
			// Distance from root
			int dist = 0;
			Node n2 = node;
			while (!n2.isRoot()) {
				
				
				// Is the other child still alive?
				boolean nodeIsReal = false;
				Node otherChild = n2.getParent().getLeft() == n2 ? n2.getParent().getRight() : n2.getParent().getLeft(); 
				for (Node leaf : otherChild.getAllLeafNodes()) {
					if (leaf.getHeight() == 0) {
						nodeIsReal = true;
						break;
					}
				}
						
				if (nodeIsReal) {
					dist ++;
				}
				
				n2 = n2.getParent();
						
			}
			leafDistancesFromRoot[node.getNr()] = dist;
			return;
		}
		
		// Children
		int inheritStateChild = random.nextFloat() < p ? stateNr : random.nextInt(n);
        sampleLabelsInClade(node.getLeft(), p, n, states, inheritStateChild, random);
        sampleLabelsInClade(node.getRight(), p, n, states, -1, random);
		
		
		
	}
	
	
	@Override
    public void init(final PrintStream out) {
        out.print(getID() + "\t" + "nstates\t");
        out.print("nstubs\t");
        for (int i = 0; i < leafDistancesFromRoot.length; i ++) {
        	String taxon = "" + i; // treeInput.get().getNode(i).getID();
        	out.print("dist" + taxon + "\t");
        	out.print("present" + taxon + "\t");
        }
        
        for (int i = 0; i < ntaxaInput.get()*2-1; i ++) {
        	out.print("nstubs.branch." + (i+1) + "\t");
        }
        
        
        
    }

    @Override
    public void log(final long sample, final PrintStream out) {
    	
    	
        out.print(getCurrentLogP() + "\t" + leafStates.size() + "\t");
        
        // Number of stubs
        out.print(getNstubs() + "\t");
        
        
        // Distance to root of each leaf, and whether it is represented in the root
        for (int i = 0; i < leafDistancesFromRoot.length; i ++) {
        	out.print(leafDistancesFromRoot[i] + "\t");
        	out.print((leafReps[i] ? 1 : 0) + "\t");
        }
        
        
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
		if (!allDescendantsAreExtinct(left) && !allDescendantsAreExtinct(right)) {
			return node;
		}
		
		
		// Left is extinct, so root is on right side
		if (allDescendantsAreExtinct(left)) {
			return getEmpiricalRoot(right);
		}
		
		// Right is extinct, so root is on left side
		else if (allDescendantsAreExtinct(right)) {
			return getEmpiricalRoot(left);
		}
		
		
		throw new IllegalArgumentException("Unexpected error: all children are extinct");
    	
    	
    }

	private int countStubs(Node node, int nstubs) {
		
		
		if (allDescendantsAreExtinct(node)) {
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


	private boolean allDescendantsAreExtinct(Node node) {
		
		
		if (node.isLeaf()) {
			if (node.getHeight() > 0) return true;
			return false;
		}
		
		
		boolean allChildrenExtinct = true;
		for (Node child : node.getChildren()) {
			if (!allDescendantsAreExtinct(child)){
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
        conditions.add(pInput.get().getID());
        conditions.add(nInput.get().getID());
        if (r0Input.get() != null) conditions.add(r0Input.get().getID());
        return conditions;
    }


	

}



