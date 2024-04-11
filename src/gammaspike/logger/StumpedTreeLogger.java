package gammaspike.logger;

import java.io.PrintStream;
import java.util.List;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.Parameter;
import beast.base.inference.parameter.RealParameter;
import gammaspike.tree.Stubs;
import beast.base.evolution.TreeWithMetaDataLogger;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;



@Description("Logs a tree and its concert model metadata")
public class StumpedTreeLogger extends TreeWithMetaDataLogger {
	
	final public Input<Stubs> stubsInput = new Input<>("stubs", "stubs model", Input.Validate.OPTIONAL);

	private boolean sortTree;
	Stubs stubs;
	
	@Override
	public void initAndValidate() {
		stubs = stubsInput.get();
		sortTree = false; // Do not sort the tree because left/right order matters 
		
		super.initAndValidate();
	}

	@Override
	public void close(PrintStream out) {
		treeInput.get().close(out);
	}

	@Override
	public void init(PrintStream out) {
		treeInput.get().init(out);
	}

	@Override
	public void log(long sample, PrintStream out) {
		
		 // make sure we get the current version of the inputs
        Tree tree = (Tree) treeInput.get().getCurrent();
        List<Function> metadata = parameterInput.get();
        for (int i = 0; i < metadata.size(); i++) {
        	if (metadata.get(i) instanceof StateNode) {
        		metadata.set(i, ((StateNode) metadata.get(i)).getCurrent());
        	}
        }
        BranchRateModel.Base branchRateModel = clockModelInput.get();
        // write out the log tree with meta data
        out.print("tree STATE_" + sample + " = ");

        if (sortTree) {
            tree.getRoot().sort();
        }

        out.print(toNewick(tree.getRoot(), metadata, branchRateModel));
        //out.print(tree.getRoot().toShortNewick(false));
        out.print(";");
		
	}
	
	

	
	
	protected String toNewick(Node node, List<Function> metadataList, BranchRateModel.Base branchRateModel) {
		
		boolean firstMetadata = true;
		StringBuffer buf = new StringBuffer();
		if (node.getLeft() != null) {
			buf.append("(");
			buf.append(toNewick(node.getLeft(), metadataList, branchRateModel));
			if (node.getRight() != null) {
				buf.append(',');
				buf.append(toNewick(node.getRight(), metadataList, branchRateModel));
			}
			buf.append(")");
		} else {
			buf.append(node.getNr() + 1);
		}
		StringBuffer buf2 = new StringBuffer();
		buf2.append("[&");
		
		// Print concertion events
		if (stubs != null) {
			int numEvents = 0;
			int totalEventCount = 0;
			List<Stubs.Stub> stubsList = stubs.getSortedStubsOnBranch(node);
			numEvents += stubsList.size();
			
			for (int eventNr = 0; eventNr < stubsList.size(); eventNr++) {
				Stubs.Stub event = stubsList.get(eventNr); 
				//System.out.println("\t\t" + event);
				event.toMetaData(buf2, eventNr+1+totalEventCount);
				buf2.append(",");
				
			}
			
			totalEventCount += stubsList.size();
				
			buf2.append(Stubs.getStubCountName() + "=");
			buf2.append(numEvents);
			firstMetadata = false;
		}
		
		
		// Branch rates
		if (branchRateModel != null) {
			if (!firstMetadata) buf2.append(",");
			buf2.append("rate=");
			appendDouble(buf2, branchRateModel.getRateForBranch(node));
			firstMetadata = false;
		}
		
		// Metadata
		if (!metadataList.isEmpty() && !firstMetadata) {
			buf2.append(",");
		}
		for (Function metadata : metadataList) {
			if (metadata instanceof Parameter<?>) {
				Parameter<?> p = (Parameter<?>) metadata;
				int dim = p.getMinorDimension1();
				if (p.getMinorDimension2() > node.getNr()) {
					buf2.append(((BEASTObject) metadata).getID());
					buf2.append('=');
					if (dim > 1) {
						buf2.append('{');
						for (int i = 0; i < dim; i++) {
							if (metadata instanceof RealParameter) {
								RealParameter rp = (RealParameter) metadata;
								appendDouble(buf2, rp.getMatrixValue(node.getNr(), i));
							} else {
								buf2.append(p.getMatrixValue(node.getNr(), i));
							}
							if (i < dim - 1) {
								buf2.append(',');
							}
						}
						buf2.append('}');
					} else {
						if (metadata instanceof RealParameter) {
							RealParameter rp = (RealParameter) metadata;
							appendDouble(buf2, rp.getArrayValue(node.getNr()));
						} else {
							buf2.append(metadata.getArrayValue(node.getNr()));
						}
					}
				} else {
				
				}
			} else {
				if (metadata.getDimension() > node.getNr()) {
					buf2.append(((BEASTObject) metadata).getID());
					buf2.append('=');
					buf2.append(metadata.getArrayValue(node.getNr()));
				}
			}
			if (buf2.length() > 2 && metadataList.indexOf(metadata) < metadataList.size() - 1) {
				buf2.append(",");
			}
		}
			
		
		buf2.append("]");
		if (buf2.length() > 3) {
			buf.append(buf2.toString());
		}
		buf.append(":");
		appendDouble(buf, node.getLength());
		return buf.toString();
	}
	
	
	private void appendDouble(StringBuffer buf, double d) {
		buf.append(d);
	}

}
