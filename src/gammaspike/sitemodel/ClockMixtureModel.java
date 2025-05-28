package gammaspike.sitemodel;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math.distribution.GammaDistribution;
import org.apache.commons.math.distribution.GammaDistributionImpl;

import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;

public class ClockMixtureModel extends SiteModel.Base {
	
	// List of clock models
	final public Input<Tree> treeInput = new Input<>("tree", "the tree (must be the same as the ones in the clock models)", Input.Validate.REQUIRED);
	final public Input<List<BranchRateModel>> branchRateModelsInput = new Input<>("branchRateModel", "a list of clock models to integrate over", new ArrayList<>());
	final public Input<RealParameter> clockWeightsInput = new Input<>("weights", "A vector of weights that each clock model receives", Input.Validate.REQUIRED);
	
	// Site rate heterogeneity
	final public Input<Integer> gammaCategoryCount = new Input<>("gammaCategoryCount", "gamma category count (default=zero for no gamma)", 0);
	final public Input<RealParameter> shapeParameterInput = new Input<>("shape", "shape parameter of gamma distribution. Ignored if gammaCategoryCount 1 or less");

	final public Input<RealParameter> muParameterInput = new Input<>("mutationRate", "mutation rate (defaults to 1.0)");
	
	int nClocks;
	int nRates;
	int nBranches;
	
	protected RealParameter muParameter;
	protected RealParameter clockWeightsParameter;
	protected RealParameter shapeParameter;
	protected List<BranchRateModel> clocks;
	
	protected boolean ratesKnown;
	
	
	// Site rates
	protected double[] siteRateHeterogeneityRates;
    protected double[] siteRateHeterogeneityProportions;
    
    // Site clocks
    protected double[] siteClockHeterogeneityBranchRates;
    protected double[] siteClockHeterogeneityProportions;
    
	
	@Override
	public void initAndValidate() {
		
		ratesKnown = false;
		
		
		// Parse site clock models
		this.clocks = branchRateModelsInput.get();
		this.clockWeightsParameter = clockWeightsInput.get();
		this.nClocks = this.clocks.size();
		if (this.nClocks == 0) {
			throw new IllegalArgumentException("Please provide at least one branchRateModel");
		}
		if (this.nClocks != this.clockWeightsParameter.getDimension()) {
			throw new IllegalArgumentException("Please ensure that weights is the same dimension as the number of branchRateModels");
		}
		this.nBranches = treeInput.get().getNodeCount();
		
		
		// Mutation rate
		this.muParameter = muParameterInput.get();
		if (muParameter == null) {
			muParameter = new RealParameter("1.0");
	    }
		
		 
		// Parse site rate models
		this.shapeParameter = shapeParameterInput.get();
		if (this.shapeParameter != null) {
			this.shapeParameter.setBounds(Math.max(this.shapeParameter.getLower(), 1.0E-3), Math.min(this.shapeParameter.getUpper(), 1.0E3));
		}
		
		
		refresh();
		
		// Do we also need to add the clock model parameters to the list??
		this.addCondition(muParameterInput);
		this.addCondition(shapeParameterInput);
		this.addCondition(clockWeightsInput);
		this.addCondition(treeInput);
		
		
	}
	
	
	
	@Override
    protected void refresh() {
		
		if (shapeParameter != null) {
            nRates = gammaCategoryCount.get();
            if (nRates < 1) {
            	if (nRates < 0) {
            		Log.warning.println("SiteModel: Invalid category count (" + nRates + ") Setting category count to 1");
            	}
                nRates = 1;
            }
        } else {
            nRates = 1;
        }


		// Site rate categories
        siteRateHeterogeneityRates = new double[nRates];
        siteRateHeterogeneityProportions = new double[nRates];
        
        // Clock categories
        siteClockHeterogeneityBranchRates = new double[nClocks * nBranches];
        siteClockHeterogeneityProportions = new double[nClocks];

        
        calculateCategoryRates();
        //ratesKnown = false;
		
		
	}
	
	
	
	public double getCategoryRate(int cat) {
		if (!ratesKnown) {
            calculateCategoryRates();
        }
		return siteRateHeterogeneityRates[cat];
	}
	
	public double getCategoryBranchRate(int cat, Node node) {
//		if (!ratesKnown) {
//            calculateCategoryRates();
//        }
//		int nodeNr = node.getNr();
//		int index = cat*this.nClocks + nodeNr;
//		return this.siteClockHeterogeneityBranchRates[index];
		
		// This is inefficient and should be cached
		BranchRateModel clockModel = this.clocks.get(cat);
		double branchRate = clockModel.getRateForBranch(node);
		return branchRate;
		
		
	}
	
//	private void setCategoryBranchRate(int cat, Node node, double value) {
//		int nodeNr = node.getNr();
//		int index = cat*this.nClocks + nodeNr;
//		this.siteClockHeterogeneityBranchRates[index] = value;
//	}

	@Override
	public boolean integrateAcrossCategories() {
		return true;
	}

	@Override
	public int getCategoryCount() {
		return nClocks*nRates;
	}

	@Override
	public int getCategoryOfSite(int site, Node node) {
		throw new IllegalArgumentException("Integrating across categories");
	}

	@Override
	public double getRateForCategory(int category, Node node) {
		
		
		int siteCategory = category % this.nRates;
		int clockCategory = (int) Math.floor(category / this.nRates); 
		
		// The rate of this category at this node is equal to the clock rate * branch rate * mu
        return getCategoryRate(siteCategory) * getCategoryBranchRate(clockCategory, node) * muParameter.getValue();
	}

	/**
     * return category rates
     *
     * @param node rates to which the rates apply. Typically, the rates will be uniform
     *             throughout the tree and the node argument is ignored.
     */
    @Override
    public double[] getCategoryRates(final Node node) {
    	final double[] rates = new double[this.getCategoryCount()];
    	for (int category = 0; category < this.getCategoryCount(); category++) {
			rates[category] = getRateForCategory(category, node);
    	}
        return rates;
    }
    
    

    

	@Override
	public double getProportionForCategory(int category, Node node) {
		
		synchronized (this) {
            if (!ratesKnown) {
                calculateCategoryRates();
            }
        }
		
		int siteCategory = category % this.nRates;
		int clockCategory = (int) Math.floor(category / this.nRates); 
		double pSite = siteRateHeterogeneityProportions[siteCategory];
		double pClock = this.clockWeightsParameter.getArrayValue(clockCategory);
		return pSite*pClock;
		
	}

	@Override
	public double[] getCategoryProportions(Node node) {
		final double[] props = new double[this.getCategoryCount()];
		for (int category = 0; category < this.getCategoryCount(); category++) {
			props[category] = getProportionForCategory(category, null);
		}
		return props;
	}
	
	
    @Override
    public void store() {
        super.store();
    } 

    @Override
    public void restore() {
        super.restore();
        ratesKnown = false;
    }
    
    
    
    @Override
    protected boolean requiresRecalculation() {
    	
        // Check whether any of the non-substitution model parameters changed
        if (this.nRates > 1 && InputUtil.isDirty(shapeParameterInput)) {
        	ratesKnown = false;
        }
        if (InputUtil.isDirty(muParameterInput)) {
        	ratesKnown = false;
        }
        if (InputUtil.isDirty(branchRateModelsInput)) {
        	ratesKnown = false;
        }	
       // if (InputUtil.isDirty(treeInput)) {
        	//ratesKnown = false;
        //}
        
        
        return true;
    }
	
	
	/**
     * discretisation of gamma distribution with equal proportions in each
     * category
     */
    protected void calculateCategoryRates() {
        double propVariable = 1.0;
        int cat = 0;


        if (shapeParameter != null) {

            final double a = shapeParameter.getValue();
            double mean = 0.0;
            final int gammaCatCount = nRates - cat;

            final GammaDistribution g = new GammaDistributionImpl(a, 1.0 / a);
            for (int i = 0; i < gammaCatCount; i++) {
                try {
                	
                    // Using beast1 style gamma
                    siteRateHeterogeneityRates[i + cat] = GammaDistributionQuantile((2.0 * i + 1.0) / (2.0 * gammaCatCount), a, 1.0 / a);

                } catch (Exception e) {
                    e.printStackTrace();
                    Log.err.println("Something went wrong with the gamma distribution calculation");
                    System.exit(-1);
                }
                mean += siteRateHeterogeneityRates[i + cat];

                siteRateHeterogeneityProportions[i + cat] = propVariable / gammaCatCount;
            }

            mean = (propVariable * mean) / gammaCatCount;

            for (int i = 0; i < gammaCatCount; i++) {

                siteRateHeterogeneityRates[i + cat] /= mean;
            }
        } else {
            siteRateHeterogeneityRates[cat] = 1.0 / propVariable;
            siteRateHeterogeneityProportions[cat] = propVariable;
        }


        ratesKnown = true;
    }
    
    
    protected double GammaDistributionQuantile(double y, double shape, double scale) {
        return 0.5 * scale * pointChi2(y, 2.0 * shape);
    }

    double pointChi2(double prob, double v) {
        // Returns z so that Prob{x<z}=prob where x is Chi2 distributed with df
        // = v
        // RATNEST FORTRAN by
        // Best DJ & Roberts DE (1975) The percentage points of the
        // Chi2 distribution. Applied Statistics 24: 385-388. (AS91)

        double e = 0.5e-6, aa = 0.6931471805, g;
        double xx, c, ch, a, q, p1, p2, t, x, b, s1, s2, s3, s4, s5, s6;

        if (prob < 0.000002 || prob > 0.999998 || v <= 0) {
            throw new IllegalArgumentException("Error SiteModel 102: Arguments out of range");
        }
        g = GammaFunctionlnGamma(v / 2);
        xx = v / 2;

        c = xx - 1;
        if (v < -1.24 * Math.log(prob)) {
            ch = Math.pow((prob * xx * Math.exp(g + xx * aa)), 1 / xx);
            if (ch - e < 0) {
                return ch;
            }
        } else {
            if (v > 0.32) {
                x = NormalDistributionQuantile(prob, 0, 1);
                p1 = 0.222222 / v;
                ch = v * Math.pow((x * Math.sqrt(p1) + 1 - p1), 3.0);
                if (ch > 2.2 * v + 6) {
                    ch = -2 * (Math.log(1 - prob) - c * Math.log(.5 * ch) + g);
                }
            } else {
                ch = 0.4;
                a = Math.log(1 - prob);

                do {
                    q = ch;
                    p1 = 1 + ch * (4.67 + ch);
                    p2 = ch * (6.73 + ch * (6.66 + ch));
                    t = -0.5 + (4.67 + 2 * ch) / p1
                            - (6.73 + ch * (13.32 + 3 * ch)) / p2;
                    ch -= (1 - Math.exp(a + g + .5 * ch + c * aa) * p2 / p1)
                            / t;
                } while (Math.abs(q / ch - 1) - .01 > 0);
            }
        }
        do {
            q = ch;
            p1 = 0.5 * ch;
            if ((t = GammaFunctionincompleteGammaP(xx, p1, g)) < 0) {
                throw new IllegalArgumentException("Error SiteModel 101: Arguments out of range: t < 0");
            }
            p2 = prob - t;
            t = p2 * Math.exp(xx * aa + g + p1 - c * Math.log(ch));
            b = t / ch;
            a = 0.5 * t - b * c;

            s1 = (210 + a * (140 + a * (105 + a * (84 + a * (70 + 60 * a))))) / 420;
            s2 = (420 + a * (735 + a * (966 + a * (1141 + 1278 * a)))) / 2520;
            s3 = (210 + a * (462 + a * (707 + 932 * a))) / 2520;
            s4 = (252 + a * (672 + 1182 * a) + c * (294 + a * (889 + 1740 * a))) / 5040;
            s5 = (84 + 264 * a + c * (175 + 606 * a)) / 2520;
            s6 = (120 + c * (346 + 127 * c)) / 5040;
            ch += t
                    * (1 + 0.5 * t * s1 - b
                    * c
                    * (s1 - b
                    * (s2 - b
                    * (s3 - b
                    * (s4 - b * (s5 - b * s6))))));
        } while (Math.abs(q / ch - 1) > e);

        return (ch);
    }
    
    

    /**
     * log Gamma function: ln(gamma(alpha)) for alpha>0, accurate to 10 decimal places
     *
     * @param alpha argument
     * @return the log of the gamma function of the given alpha
     */
    double GammaFunctionlnGamma(final double alpha) {
        // Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
        // Communications of the Association for Computing Machinery, 9:684

        double x = alpha, f = 0.0, z;

        if (x < 7) {
            f = 1;
            z = x - 1;
            while (++z < 7) {
                f *= z;
            }
            x = z;
            f = -Math.log(f);
        }
        z = 1 / (x * x);

        return
                f + (x - 0.5) * Math.log(x) - x + 0.918938533204673 +
                        (((-0.000595238095238 * z + 0.000793650793651) *
                                z - 0.002777777777778) * z + 0.083333333333333) / x;
    }
	
    

    /**
     * Incomplete Gamma function P(a,x) = 1-Q(a,x)
     * (a cleanroom implementation of Numerical Recipes gammp(a,x);
     * in Mathematica this function is 1-GammaRegularized)
     *
     * @param a        parameter
     * @param x        argument
     * @param lnGammaA precomputed lnGamma(a)
     * @return function value
     */
    double GammaFunctionincompleteGammaP(double a, double x, double lnGammaA) {
        return incompleteGamma(x, a, lnGammaA);
    }


    /**
     * Returns the incomplete gamma ratio I(x,alpha) where x is the upper
     * limit of the integration and alpha is the shape parameter.
     *
     * @param x              upper limit of integration
     * @param alpha          shape parameter
     * @param ln_gamma_alpha the log gamma function for alpha
     * @return the incomplete gamma ratio
     */
    double incompleteGamma(double x, double alpha, double ln_gamma_alpha) {
        // (1) series expansion     if (alpha>x || x<=1)
        // (2) continued fraction   otherwise
        // RATNEST FORTRAN by
        // Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
        // 19: 285-287 (AS32)

        double accurate = 1e-8, overflow = 1e30;
        double factor, gin, rn, a, b, an, dif, term;
        double pn0, pn1, pn2, pn3, pn4, pn5;

        if (x == 0.0) {
            return 0.0;
        }
        if (x < 0.0 || alpha <= 0.0) {
            throw new IllegalArgumentException("Error SiteModel 103: Arguments out of bounds");
        }

        factor = Math.exp(alpha * Math.log(x) - x - ln_gamma_alpha);

        if (x > 1 && x >= alpha) {
            // continued fraction
            a = 1 - alpha;
            b = a + x + 1;
            term = 0;
            pn0 = 1;
            pn1 = x;
            pn2 = x + 1;
            pn3 = x * b;
            gin = pn2 / pn3;

            do {
                a++;
                b += 2;
                term++;
                an = a * term;
                pn4 = b * pn2 - an * pn0;
                pn5 = b * pn3 - an * pn1;

                if (pn5 != 0) {
                    rn = pn4 / pn5;
                    dif = Math.abs(gin - rn);
                    if (dif <= accurate) {
                        if (dif <= accurate * rn) {
                            break;
                        }
                    }

                    gin = rn;
                }
                pn0 = pn2;
                pn1 = pn3;
                pn2 = pn4;
                pn3 = pn5;
                if (Math.abs(pn4) >= overflow) {
                    pn0 /= overflow;
                    pn1 /= overflow;
                    pn2 /= overflow;
                    pn3 /= overflow;
                }
            } while (true);
            gin = 1 - factor * gin;
        } else {
            // series expansion
            gin = 1;
            term = 1;
            rn = alpha;
            do {
                rn++;
                term *= x / rn;
                gin += term;
            }
            while (term > accurate);
            gin *= factor / alpha;
        }
        return gin;
    }

    double NormalDistributionQuantile(double z, double m, double sd) {
        return m + Math.sqrt(2.0) * sd * ErrorFunctionInverseErf(2.0 * z - 1.0);
    }

    /**
     * inverse error function
     *
     * @param z argument
     * @return function value
     */
    double ErrorFunctionInverseErf(double z) {
        return ErrorFunctionPointNormal(0.5 * z + 0.5) / Math.sqrt(2.0);
    }


    // Private

    // Returns z so that Prob{x<z}=prob where x ~ N(0,1) and (1e-12) < prob<1-(1e-12)

    double ErrorFunctionPointNormal(double prob) {
        // Odeh RE & Evans JO (1974) The percentage points of the normal distribution.
        // Applied Statistics 22: 96-97 (AS70)

        // Newer methods:
        // Wichura MJ (1988) Algorithm AS 241: the percentage points of the
        // normal distribution.  37: 477-484.
        // Beasley JD & Springer SG  (1977).  Algorithm AS 111: the percentage
        // points of the normal distribution.  26: 118-121.

        double a0 = -0.322232431088, a1 = -1, a2 = -0.342242088547, a3 = -0.0204231210245;
        double a4 = -0.453642210148e-4, b0 = 0.0993484626060, b1 = 0.588581570495;
        double b2 = 0.531103462366, b3 = 0.103537752850, b4 = 0.0038560700634;
        double y, z, p1;

        p1 = (prob < 0.5 ? prob : 1 - prob);
        if (p1 < 1e-20) {
            throw new IllegalArgumentException("Error SiteModel 104: Argument prob out of range");
        }

        y = Math.sqrt(Math.log(1 / (p1 * p1)));
        z = y + ((((y * a4 + a3) * y + a2) * y + a1) * y + a0) / ((((y * b4 + b3) * y + b2) * y + b1) * y + b0);
        return (prob < 0.5 ? -z : z);
    }
	

}
