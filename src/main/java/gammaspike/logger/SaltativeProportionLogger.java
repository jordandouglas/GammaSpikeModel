package gammaspike.logger;

import java.io.PrintStream;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import gammaspike.clockmodel.PunctuatedRelaxedClockModel;

@Description("Logs the proportion of evolutionary distance across the tree that can be explained by abrupt change, as opposed to gradual")
public class SaltativeProportionLogger  extends BEASTObject implements Loggable {

	final public Input<Tree> treeInput = new Input<>("tree", "tree", Validate.REQUIRED);
	final public Input<PunctuatedRelaxedClockModel> clockInput = new Input<>("clock", "the abrupt clock model", Validate.REQUIRED);

	@Override
	public void initAndValidate() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void init(PrintStream out) {
		out.print(this.clockInput.get().getID() + ".ProportionOfSaltation\t");
		
	}

	@Override
	public void log(long sample, PrintStream out) {
		
		
		double totalChange = 0;
		double abruptChange = 0;
		
		for (Node node : treeInput.get().getNodesAsArray()) {
			
			if (node.getLength() <= 0 || node.isDirectAncestor() || node.isRoot()) continue;
			
			// Total distance on branch
			double branchTime = node.isRoot() ? 0 : node.getLength();
			double branchRate = clockInput.get().getRateForBranch(node);
			double branchDistance = branchTime*branchRate;
			
			
			totalChange += branchDistance;
			
			
			// Abrupt change
			double abruptChangeBranch = clockInput.get().getBurstSize(node);
			abruptChange += abruptChangeBranch;
			
			
		}
		
		out.print(abruptChange/totalChange + "\t");
		
	}

	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
		
	}

}
