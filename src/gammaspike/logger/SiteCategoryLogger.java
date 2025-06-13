package gammaspike.logger;

import java.io.PrintStream;

import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.core.Loggable;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.likelihood.LikelihoodCore;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.sitemodel.SiteModelInterface.Base;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;

public class SiteCategoryLogger extends BEASTObject implements Loggable {
	
	final public Input<TreeLikelihood> likelihoodInput = new Input<>("likelihood", "Tree likelihood", Input.Validate.REQUIRED);
	
	int siteCount;
	int nrOfPatterns;
	int categoryCount;
	
	@Override
	public void initAndValidate() {
		
		Alignment data = likelihoodInput.get().dataInput.get();
		this.siteCount = data.getSiteCount();
		this.nrOfPatterns = data.getPatternCount();
		
		SiteModel.Base sm = (Base) likelihoodInput.get().siteModelInput.get();
		this.categoryCount = sm.getCategoryCount();
		
	}

	@Override
	public void init(PrintStream out) {
		
		for (int i = 0; i < siteCount; i++) {
			out.print((getID() != null ? getID() : "site")  + "." + (i+1) + "\t");  
		}
		
	}

	@Override
	public void log(long sample, PrintStream out) {
		
		
		// Likelihood
		TreeLikelihood likelihood = likelihoodInput.get();
		Tree tree = (Tree) likelihood.treeInput.get();
		//tree.setEverythingDirty(true);
		double logP = likelihood.calculateLogP();
		
		// Get category weights
		SiteModel.Base sm = (Base) likelihood.siteModelInput.get();
		double[] categoryProportions = sm.getCategoryProportions(null);
		
		// Get the partials at the root for each category
		double[] rootPartials = new double[this.nrOfPatterns*this.categoryCount];
		LikelihoodCore core = likelihood.getLikelihoodCore();
		core.getNodePartials(tree.getRoot().getNr(), rootPartials);
		
		
		
		
		// Posterior vector
		double[][] posteriorOfEachCategoryPerPattern = new double[this.nrOfPatterns][this.categoryCount];
		double[] probsPattern = new double[this.categoryCount];
		int k = 0;
		for (int patternNum = 0; patternNum < this.nrOfPatterns; patternNum ++) {
			
			
			// Posterior probs
			double pSum = 0;
			for (int categoryNum = 0; categoryNum < this.categoryCount; categoryNum++) {
				//Log.warning("patternNum=" + patternNum + " , categoryNum=" + categoryNum + " " + rootPartials[k] + "  / " + categoryProportions[categoryNum]);
				probsPattern[categoryNum] = rootPartials[k] * categoryProportions[categoryNum];
				pSum += probsPattern[categoryNum];
				k++;
			}
			
			// Normalise
			for (int categoryNum = 0; categoryNum < this.categoryCount; categoryNum++) {
				posteriorOfEachCategoryPerPattern[patternNum][categoryNum] = probsPattern[categoryNum] / pSum;
			}
			
			
			
		}
        
		// We don't need these
//		double[] rootFrequencies = likelihood.getSubstitutionModel().getFrequencies();
//        if (likelihood.rootFrequenciesInput.get() != null) {
//        	rootFrequencies = likelihood.rootFrequenciesInput.get().getFreqs();
//        }
		
		
		
		for (int i = 0; i < siteCount; i++) {
			
			int patternNum = likelihood.dataInput.get().getPatternIndex(i);
			double[] pvector = posteriorOfEachCategoryPerPattern[patternNum];
			double sum = 0;
			for (int j = 0 ; j < pvector.length; j ++) {
				sum += pvector[j];
			}
			if (sum < 1 - 1e-6 || sum > 1 + 1e-6) { // Should equal 1
				out.print(-1 + "\t");  
			}else {
				out.print(Randomizer.randomChoicePDF(pvector) + "\t");  
			}
			
			
		}
		
	}

	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
		
	}
	
	

}



