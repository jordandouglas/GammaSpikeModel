package gammaspike.flabel;

import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.CalculationNode;

public class FlabelData extends CalculationNode {
	
	final public Input<TreeInterface> treeInput = new Input<>("tree", "the tree", Input.Validate.REQUIRED);
	final public Input<String> taxonInput = new Input<>("taxon", "the taxon in the tree", Input.Validate.REQUIRED);
	final public Input<Integer> labelInput = new Input<>("val", "the value of this", Input.Validate.REQUIRED);
	

	@Override
	public void initAndValidate() {
		
		
		// Check taxon exists
		boolean foundTaxon = false;
		for (Node leaf : treeInput.get().getNodesAsArray()) {
			if (!leaf.isLeaf()) continue;
			//Log.warning("LEAF: '" + leaf.getID() + "'|'" + taxonInput.get() + "'");
			if (leaf.getID().equals(taxonInput.get())) {
				
				foundTaxon = true;
				break;
			}
		}
		
		
		if (!foundTaxon) {
			throw new IllegalArgumentException("Cannot find " + taxonInput.get() + " in the tree");
		}
		
		
	}
	
	public String getTaxon() {
		return taxonInput.get();
	}
	
	public int getLabel() {
		return labelInput.get();
	}

}
