package gammaspike.util;


import java.io.PrintStream;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Runnable;
import beastfx.app.tools.Application;
import beastfx.app.treeannotator.TreeAnnotator;
import beastfx.app.treeannotator.TreeAnnotator.MemoryFriendlyTreeSet;
import beastfx.app.util.OutFile;
import beastfx.app.util.TreeFile;

@Description("Extract effective rates from tree log produced by gamma spike model")
public class GetEffectiveRates extends Runnable {
	final public Input<TreeFile> treesInput = new Input<>("trees","NEXUS file containing a tree set", Validate.REQUIRED);
	final public Input<Integer> burnInPercentageInput = new Input<>("burnin", "percentage of trees to used as burn-in (and will be ignored)", 10);
	final public Input<OutFile> outputInput = new Input<>("out","output file. Print to stdout if not specified");

	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		// get trees
		MemoryFriendlyTreeSet srcTreeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(treesInput.get().getPath(), burnInPercentageInput.get());
		srcTreeSet.reset();

		PrintStream out = System.out;
        if (outputInput.get() != null) {
			Log.warning("Writing to file " + outputInput.get().getPath());
        	out = new PrintStream(outputInput.get());
        }
        
        // print header
		Tree tree = srcTreeSet.next();
		out.print("Sample\t");
		for (int i = 0; i < tree.getLeafNodeCount(); i++) {
			out.print("eRate" + (i+1) + "\t");
		}
		out.println();
		
		// process trees
        srcTreeSet.reset();
        int k = 1;
        while (srcTreeSet.hasNext()) {
        	tree = srcTreeSet.next();
        	Node [] nodes = tree.getNodesAsArray();
        	out.print(k+"\t");
    		for (int i = 0; i < tree.getLeafNodeCount(); i++) {
    			double rate = meta(nodes[i], "branchRates");
    			double spike = meta(nodes[i], "spikes");
    			double len = nodes[i].getLength();
    			out.print((rate + spike/len) + "\t");
    		}
    		out.println();
    		k++;
		}
		
        //close
        if (outputInput.get() != null) {
        	out.close();
        }
        Log.warning.println("Done.");
	}

	private double meta(Node node, String str) {
		Object o = (Double) node.getMetaData(str);
		if (o instanceof String) {
			double value = Double.valueOf((String) o);
			return value;
		}
		if (o instanceof Double) {
			return (Double) o;
		}
		return 0;
	}

	public static void main(String[] args) throws Exception {
		new Application(new GetEffectiveRates(), "Get effectiver rates from spike model tree log", args);

	}

}
