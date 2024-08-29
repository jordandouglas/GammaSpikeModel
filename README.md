# GammaSpikeModel

The gamma spike model is a phylogenetic clock model for BEAST 2. This method tests for punctuated equilibrium, by assigning each branch a gradual clock-rate and an instaneaneous spike of abrupt evolution. The magnitudes of each branch rate and branch spike are independent and identically distributed.  

The genetic distance of a branch is calculated as:
```
distance = rate x length + spike
```


Advantages over other clock models:

1. Accounts for the common scenario where evolution is accelerated at the time of branching.
2. The hypothesis of punctuated equilibrium is tested during MCMC using model averaging.
3. If the hypothesis is rejected, then MCMC falls back on the relaxed clock ([ORC](https://github.com/jordandouglas/ORC) model).
4. If the hypothesis is accepted, then the estimated tree is likely to be more accurate.
5. In either case, the number of unobserved speciation events on each branch is also estimated (i.e., the *stubs*).



## Installation instructions


This package requires BEAST 2.7.7. or newer.

1. Launch BEAUti
2. Click on `File` -> `Manage Packages`
3. Install GammaSpikeModel. If it is not in the list of packages, you may need to add an extra package repository as follows:
- Click the `Packager repositories` button. A dialog pops up.
- Click the Add URL button. A dialog is shown where you can enter https://raw.githubusercontent.com/CompEvol/CBAN/master/packages-extra-2.7.xml
- Click `OK`. There should be an extra entry in the list.
- Click `Done`.
- After a short delay, GammaSpikeModel should appear in the list of packages.



## Setting up an analysis


1. Launch BEAUti.


2. Import data and set up site model, as per usual.

3. Open the `Clock Model` tab and select `Punctuated Relaxed Clock Model`.
	- This will introduce the following parameters
		- The mean spike size `spikeMean` (under a Gamma distribution). High mean = larger expected spike sizes. 
		- The shape of spike sizes `spikeShape` (Gamma distribution). High shape = smaller variance of spike sizes.
		- One `spike` for every branch in the tree, whose sizes are Gamma distributed under the prior.
		- A boolean model indicator `useSpikeModel` that determines whether the spikes are being used or not.  

4. Open the `Priors` tab and select the `Stumped Tree Prior`. 
	- This uses the fossilised birth-death tree prior with the following parameters
		- Birth rate `lambdaXX` .
		- Reproduction number `R0`, which is assumed to be greater than 1. 
		- Sampling proportion `XXXX`. If all taxa are extant, then set this to 0 and the model is just a birth-death process.
	- By using this prior, the number of stubs on each branch will be logged and inform the clock model spike sizes.
	- Sampled ancestors are estimated. 
	- If a stumped tree prior is not selected, the clock model will assume there are no stubs on any branch.
	- Currently, stubs are only available for the fossilised-birth-death model (and not coalescent or skyline models).

5. Configure priors and save the XML file as per usual.




## Spikes
The total sum of spikes is estimated for each branch. The total spike sum `Si` on branch `i` follows a Gamma distribution under the prior:

```
Si ~ Gamma(shape=spikeShape*(ni + 1), scale=spikeMean/spikeShape)
```

where `ni` is the number of stubs on branch `i`. Longer branches and older branches usually have more stubs. The reported `spike` on each branch is the total sum of all `ni + 1` spikes along that branch.

Then `spikeMean` is the average rate that sites change at each speciation event (observed or unobserved). For example an average spike size of 0.01 means that 1% of all sites are expected to change (possibly back into the original state) at each bifurcation. We examined 9 empirical datasets with support for punctuated equilibrium, and found that `spikeMean` estimates ranged from 0.001 to 0.07. The default prior for `spikeMean` is centered around this interval. This is an important parameter, and testing for sensitivity is recommended.



## Estimating stubs



## Hypothesis testing

The `useSpikeModel` parameter can be used for hypothesis testing. When this parameter is 1, the spike model is being used, and when 0 the relaxed clock. Every dataset is different, with some strongly favouring one model over the other, and others being uncertain. If the average value of `useSpikeModel` is over 0.9, there is strong support in favour of punctuated equilibrium. 


## Convergence during MCMC

In most instances, we found this model converges quite well despite its large parameter space  However, the `useSpikeModel` parameter can lead to a bimodal distribution that causes mixing issues.  If this model is too slow to converge, we recommend the following options:

1. Try using [Coupled MCMC](https://www.beast2.org/2020/01/14/metropolis-coupled-mcmcmc3-works.html).

2. Try running separate analyses with `useSpikeModel` respectivley fixed at either 0 or 1. Unfortunately, this configuration does not enable hypothesis testing (Bayesian model averaging). 


## Support

BEAST user forums [https://groups.google.com/g/beast-users](https://groups.google.com/g/beast-users)

Jordan Douglas jordan.douglas@auckland.ac.nz


## References

