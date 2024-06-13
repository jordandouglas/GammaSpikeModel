package gammaspike.logger;

import java.io.PrintStream;

import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Tree;

public class TaxonCountLogger extends BEASTObject implements Loggable {

	final public Input<Tree> treeInput = new Input<>("tree", "tree to be logged", Validate.REQUIRED);
	
	@Override
	public void init(PrintStream out) {
		out.print(treeInput.get().getID() + ".ntaxa\t");
	}

	@Override
	public void log(long sample, PrintStream out) {
		out.print(treeInput.get().getLeafNodeCount() + "\t");
	}

	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void initAndValidate() {
		// TODO Auto-generated method stub
		
	}

}
