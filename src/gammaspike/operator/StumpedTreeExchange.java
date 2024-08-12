package gammaspike.operator;

import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.operator.Exchange;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.StateNode;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import gammaspike.tree.Stubs;
import sa.evolution.operators.SAExchange;


@Description("Exchange operator that accounts for stub placement")
public class StumpedTreeExchange extends SAExchange {
	
	final public Input<Stubs> stubsInput = new Input<>("stubs", "stubs for the tree", Input.Validate.REQUIRED);
	
	
	@Override
    public double proposal() {
		
		// Cache branch lengths before making proposal
		Stubs stubs = stubsInput.get();
        double[] cachedBranchLengths = stubs.prepareJacobian();
		
		double HR = super.proposal();
		
		// Jacobian. Relative stub heights stay the same but absolute heights change
        double logJacobian = stubs.getLogJacobian(cachedBranchLengths);
        
        return HR + logJacobian;
		
	}

	
	
//
//
//	@Override
//	public void initAndValidate() {
//		// TODO Auto-generated method stub
//		
//	}
//
//	@Override
//	public double proposal() {
//		
//		final Tree tree = (Tree) InputUtil.get(treeInput, this);
//
//        double logHastingsRatio = 0;
//
//        if (isNarrowInput.get()) {
//            logHastingsRatio = narrow(tree);
//        } else {
//            logHastingsRatio = wide(tree);
//        }
//
//        return logHastingsRatio;
//	}
//	
//      private int isg(final Node n) {
//        return (n.getLeft().isLeaf() && n.getRight().isLeaf()) ? 0 : 1;
//      }
//
//      private int sisg(final Node n) {
//          return n.isLeaf() ? 0 : isg(n);
//      }
//
//    /**
//     * WARNING: Assumes strictly bifurcating beast.tree.
//     */
//    public double narrow(final Tree tree) {
//
//    	
//    	Stubs stubs = stubsInput.get();
//    	
//    	
//    	// Same as narrow exchange
//        final int internalNodes = tree.getInternalNodeCount();
//        if (internalNodes <= 1) {
//            return Double.NEGATIVE_INFINITY;
//        }
//        
//        if (internalNodes == 1 && tree.getRoot().isFake()) {
//            return Double.NEGATIVE_INFINITY;
//        }
//
//        Node grandParent = tree.getNode(internalNodes + 1 + Randomizer.nextInt(internalNodes));
//        while ((grandParent.getLeft().isLeaf() && grandParent.getRight().isLeaf()) || grandParent.isFake()) {
//            grandParent = tree.getNode(internalNodes + 1 + Randomizer.nextInt(internalNodes));
//        }
//
//        Node parentIndex = grandParent.getLeft();
//        Node uncle = grandParent.getRight();
//        if (parentIndex.getHeight() < uncle.getHeight()) {
//            parentIndex = grandParent.getRight();
//            uncle = grandParent.getLeft();
//        }
//
//        // Tree with dated tips
//        if (parentIndex.isLeaf() ) {
//            return Double.NEGATIVE_INFINITY;
//        }
//        
//        // Cache branch lengths before making proposal
//        double[] cachedBranchLengths = stubs.prepareJacobian();
//
//        int validGP = 0;
//        {
//            for(int i = internalNodes + 1; i < 1 + 2*internalNodes; ++i) {
//                validGP += isg(tree.getNode(i));
//            }
//        }
//
//        final int c2 = sisg(parentIndex) + sisg(uncle);
//
//        final Node child = Randomizer.nextBoolean() ? parentIndex.getLeft() : parentIndex.getRight();
//        
//        double childL1 = child.getLength();
//        double uncleL1 = uncle.getLength();
//        exchangeNodes(child, uncle, parentIndex, grandParent);
//        double childL2 = child.getLength();
//        double uncleL2 = uncle.getLength();
//
//        final int validGPafter = validGP - c2 + sisg(parentIndex) + sisg(uncle);
//
//        
//        // Now we need to adjust the stubs on the two branches that moved
//        //Log.warning("child length went from " + childL1 + " to " + childL2);
//        //Log.warning("uncle length went from " + uncleL1 + " to " + uncleL2);
//        
//        
//        double stubJacobian = stubs.getLogJacobian(cachedBranchLengths);
//        /*
//        boolean swapStubs = false;
//        double stubJacobian = 0;
//        if (swapStubs) {
//        	
//        	for (int stubNr = 0; stubNr < stubs.getStubDimension(); stubNr++) {
//            	if (!stubs.includeStub(stubNr)) continue;
//            	int stubBranchNr = stubs.getBranches().getValue(stubNr);
//            	
//            	// Swap stub indices
//            	if (stubBranchNr == child.getNr()) {
//            		stubs.getBranches().setValue(stubNr, uncle.getNr());
//            	}else if (stubBranchNr == uncle.getNr()) {
//            		stubs.getBranches().setValue(stubNr, child.getNr());
//            	}
//        	
//            }
//        	
//        }else {
//        	
//        	  
//            //RealParameter heights = stubs.getStubHeights();
//            for (int stubNr = 0; stubNr < stubs.getStubDimension(); stubNr++) {
//            	if (!stubs.includeStub(stubNr)) continue;
//            	int stubBranchNr = stubs.getBranches().getValue(stubNr);
//            	
//            	// Since relative times stay constant, but absolute change, we adjust the Jacobians
//            	if (stubBranchNr == child.getNr()) {
//            		stubJacobian += Math.log(childL2) - Math.log(childL1);
//            	}else if (stubBranchNr == uncle.getNr()) {
//            		stubJacobian += Math.log(uncleL2) - Math.log(uncleL1);
//            	}
//        	
//            }
//        
//        
//      
//        }
//        */
//        
//        
//        
//        return stubJacobian + Math.log((float)validGP/validGPafter);
//    }
//    
//    
//
//    /**
//     * WARNING: Assumes strictly bifurcating beast.tree.
//     * @param tree
//     */
//    public double wide(final Tree tree) {
//
//        final int nodeCount = tree.getNodeCount();
//        Stubs stubs = stubsInput.get();
//        
//        if (nodeCount == 3 && tree.getRoot().isFake()) {
//            return Double.NEGATIVE_INFINITY;
//        }
//        
//        // Cache branch lengths before making proposal
//        double[] cachedBranchLengths = stubs.prepareJacobian();
//
//        Node i = tree.getRoot();
//
//        while (i.isRoot() || i.isDirectAncestor()) {
//            i = tree.getNode(Randomizer.nextInt(nodeCount));
//        }
//
//        Node j = i;
//        while (j.getNr() == i.getNr() || j.isRoot() || j.isDirectAncestor()) {
//            j = tree.getNode(Randomizer.nextInt(nodeCount));
//        }
//
//        final Node p = i.getParent();
//        final Node jP = j.getParent();
//
//        if ((p != jP) && (i != jP) && (j != p) && (j.getHeight() < p.getHeight()) && (i.getHeight() < jP.getHeight())) {
//        	
//			 double iL1 = i.getLength();
//			 double jL1 = j.getLength();
//			
//			exchangeNodes(i, j, p, jP);
//			
//			 double iL2 = i.getLength();
//			 double jL2 = j.getLength();
//
//            // All the nodes on the path from i/j to the common ancestor of i/j parents had a topology change,
//            // so they need to be marked FILTHY.
//            if( markCladesInput.get() ) {
//                Node iup = p;
//                Node jup = jP;
//                while (iup != jup) {
//                    if( iup.getHeight() < jup.getHeight() ) {
//                        assert !iup.isRoot();
//                        iup = iup.getParent();
//                        iup.makeDirty(Tree.IS_FILTHY);
//                    } else {
//                        assert !jup.isRoot();
//                        jup = jup.getParent();
//                        jup.makeDirty(Tree.IS_FILTHY);
//                    }
//                }
//            }
//            
//            
//            
//            
//            double stubJacobian = stubs.getLogJacobian(cachedBranchLengths);
//            /*
//            double stubJacobian = 0;
//            //RealParameter heights = stubs.getStubHeights();
//            for (int stubNr = 0; stubNr < stubs.getStubDimension(); stubNr++) {
//            	if (!stubs.includeStub(stubNr)) continue;
//            	int stubBranchNr = stubs.getBranches().getValue(stubNr);
//            	
//            	// Since relative times stay constant, but absolute change, we adjust the Jacobians
//            	if (stubBranchNr == i.getNr()) {
//            		stubJacobian += Math.log(iL2) - Math.log(iL1);
//            	}else if (stubBranchNr == j.getNr()) {
//            		stubJacobian += Math.log(jL2) - Math.log(jL1);
//            	}
//            }
//            */
//            
//            return stubJacobian;
//        }
//
//        // Randomly selected nodes i and j are not valid candidates for a wide exchange.
//        // reject instead of counting (like we do for narrow).
//        return Double.NEGATIVE_INFINITY;
//    }
//
//	
//	
//	@Override
//    public List<StateNode> listStateNodes() {
//		Stubs stubs = stubsInput.get();
//        final List<StateNode> list = super.listStateNodes();
//        list.add(stubs.getBranches());
//        return list;
//    }

}
