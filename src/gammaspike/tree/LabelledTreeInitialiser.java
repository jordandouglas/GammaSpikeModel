package gammaspike.tree;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.coalescent.RandomTree;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;
import gammaspike.flabel.Flabel;


@Description("Initalises the tree such that each label is monophyletic. Not a great initial state, but at least it's legal")
public class LabelledTreeInitialiser extends Tree implements StateNodeInitialiser {

	
	final public Input<Flabel> flabelsInput = new Input<>("flabel", "leaf and internal labels", Input.Validate.REQUIRED);
	
	
	final double dt = 0.01;
	
	 @Override
	 public void initAndValidate() {
		 initStateNodes();
		 super.initAndValidate();
	 }
	
	
	
	 @Override
	 public void initStateNodes() {
		 
		
		// Make one clade per label
		Tree tree = m_initial.get();
		Flabel flabel = flabelsInput.get();
		
		
		// Get list of clades/labels
		Map<Integer, List<String>> classes = new HashMap<>();
		for (Node leaf : tree.getNodesAsArray()) {
			if (!leaf.isLeaf()) continue;
			
			
			int label = flabel.getLeafLabel(leaf);
			if (classes.containsKey(label)) {
				classes.get(label).add(leaf.getID());
			}else {
				List<String> array = new ArrayList<>();
				array.add(leaf.getID());
				classes.put(label, array);
			}
			
		}
		
		
		
		
		
		// For each label make a clade (caterpillar)
		List<Node> clades = new ArrayList<>();
		for (int label : classes.keySet()) {
			
			
			
			List<String> nodes = classes.get(label);

			// System.out.println("label" + label + ":" + Arrays.toString(taxa.toArray()));
		
			Node root = new Node();
			if (nodes.size() == 1) {
				root.setHeight(0);
				root.setID(nodes.get(0));
			}else {
				
				
				double height = dt;
				
				// Make 2 leaves with the right IDs
				String label1 = nodes.get(0);
				String label2 = nodes.get(1);
				Node n1 = new Node();
				n1.setHeight(0);
				n1.setID(label1);
				Node n2 = new Node();
				n2.setHeight(0);
				n2.setID(label2);
				
				// Make a caterpillar
				root.addChild(n1);
				root.addChild(n2);
				//n1.removeAllChildren(false);
				//n2.removeAllChildren(false);
				root.setHeight(height);
				height += dt;
				for (int i = 2; i < nodes.size(); i++) {
					
					// Make a leaf
					String labelN = nodes.get(i);
					Node node = new Node();
					node.setHeight(0);
					node.setID(labelN);
					
					
					// Attach it to caterpillar
					Node newRoot = new Node();
					newRoot.addChild(root);
					newRoot.addChild(node);
					newRoot.setHeight(height);
					height += dt;
					
					root = newRoot;
					
				}
				
			}
			
			
			// Store the root
			//Log.warning(label  + "=" + root);
			clades.add(root);
			
			
		}
		
		
		
		
		// Join the clades together in a cateprillar
		Node root = null;
		if (clades.size() == 1) {
			root = clades.get(0);
		}else {
			
			
			
			
			// Make a caterpillar
			Node n1 = clades.get(0);
			Node n2 = clades.get(1);
			
			
			root = new Node();
			root.addChild(n1);
			root.addChild(n2);
			
			double height = Math.max(n1.getHeight(), n2.getHeight()) + dt;
			root.setHeight(height);
			
			for (int i = 2; i < clades.size(); i++) {
				
				
				
				Node node = clades.get(i);
				Node newRoot = new Node();
				newRoot.addChild(root);
				newRoot.addChild(node);
				height = Math.max(root.getHeight(), node.getHeight()) + dt;
				newRoot.setHeight(height);
				root = newRoot;
				
			}
			
			
		}
		
		
		// Renumber leaves
		int nr = 0;
		for (Node leaf : root.getAllLeafNodes()) {
			leaf.setNr(nr);
			nr++;
		}
		
		// Renumber internal
		for (Node node : root.getAllChildNodesAndSelf()) {
			if (node.isRoot()) continue;
			if (node.isLeaf()) continue;
			node.setNr(nr);
			nr++;
		}
		root.setNr(nr);
				
		
		setRoot(root);
        initArrays();
        m_initial.get().assignFromWithoutID(this);
        
        flabel.initialiseIndicators(true);
        flabel.resetLeafLabels();
		
		// Make tree
		//m_initial.get().setRoot(root);
		//m_initial.get().initArrays();
		
		//Log.warning(m_initial.get().getRoot().toNewick());
		
		//Tree newTree = new Tree();
		//newTree.setRoot(root);
		//newTree.initArrays();
		//m_initial.get().assignFromWithoutID(newTree);
		
		
	}
	
	
	
	@Override
	public void getInitialisedStateNodes(List<StateNode> stateNodes) {
		stateNodes.add(m_initial.get());
	}
	
	
}


