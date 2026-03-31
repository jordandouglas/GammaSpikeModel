package gammaspike.logger;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.StateNode;
import gammaspike.tree.Stubs;
import beast.base.spec.evolution.TreeWithMetaDataLogger;
import beast.base.spec.evolution.branchratemodel.Base;
import beast.base.spec.type.Tensor;
import beast.base.spec.type.Vector;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;



@Description("Logs a tree and its concert model metadata")
public class StumpedTreeLogger extends TreeWithMetaDataLogger {
	
	final public Input<Stubs> stubsInput = new Input<>("stubs", "stubs model", Input.Validate.OPTIONAL);
	final public Input<Boolean> printStubLocationsInput = new Input<>("printStubLocations", "should stub heights be printed too", false);

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
        List<Tensor<?,?>> metadata = parameterInput.get();
        List<Tensor<?,?>> currentmetadata = new ArrayList<>();

        for (int i = 0; i < metadata.size(); i++) {
        	if (metadata.get(i) instanceof StateNode sn) {
        		currentmetadata.add((Tensor<?,?>) sn.getCurrent());
        	} else {
        		currentmetadata.add(metadata.get(i));
        	}
        }
        Base branchRateModel = clockModelInput.get();
        // write out the log tree with meta data
        out.print("tree STATE_" + sample + " = ");

        if (sortTree) {
            tree.getRoot().sort();
        }

        out.print(toNewick(tree.getRoot(), currentmetadata, branchRateModel, sample));
        out.print(";");
    }
	

	
	
	protected String toNewick(Node node, List<Tensor<?, ?>> metadataList, Base branchRateModel, long sampleNr) {
		
		boolean firstMetadata = true;
		StringBuffer buf = new StringBuffer();
		if (node.getLeft() != null) {
			buf.append("(");
			buf.append(toNewick(node.getLeft(), metadataList, branchRateModel, sampleNr));
			if (node.getRight() != null) {
				buf.append(',');
				buf.append(toNewick(node.getRight(), metadataList, branchRateModel, sampleNr));
			}
			buf.append(")");
		} else {
			buf.append(node.getNr() + 1);
		}
		StringBuffer buf2 = new StringBuffer();
		buf2.append("[&");
		
		// Print stubs
		if (stubs != null && !node.isRoot()) {
			int numEvents = 0;
			int totalEventCount = 0;
			
			if (!stubs.estimateStubs()) {
				int nstubs = stubs.sampleNStubsOnBranch(node.getNr(), sampleNr);
				numEvents += nstubs;
				totalEventCount += nstubs;
			}else if (!stubs.getReversibleJump()) {
				int nstubs = stubs.getNStubsOnBranch(node.getNr());
				numEvents += nstubs;
				totalEventCount += nstubs;
			}else {
			
			
				List<Stubs.Stub> stubsList = stubs.getSortedStubsOnBranch(node);
				numEvents += stubsList.size();
				
				if (printStubLocationsInput.get()) {
					for (int eventNr = 0; eventNr < stubsList.size(); eventNr++) {
						Stubs.Stub event = stubsList.get(eventNr); 
						//System.out.println("\t\t" + event);
						event.toMetaData(buf2, eventNr+1+totalEventCount);
						buf2.append(",");
						
					}
				}
				
				totalEventCount += stubsList.size();
			
			}
				
			//buf2.append(Stubs.getStubCountName() + "=");
			buf2.append("nstubs=");
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
		
		boolean needsComma = false;
		for (Tensor<?,?> metadata : metadataList) {
			if (metadata instanceof Vector) {
				
				// TODO hadle matrices/vectors
				needsComma = true;
			} else {
				if (metadata.size() > node.getNr()) {
					if (needsComma) {
						buf2.append(",");
					}
					buf2.append(((BEASTObject) metadata).getID());
					buf2.append('=');
					buf2.append(metadata.get(node.getNr()));
					needsComma = true;
				}
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
