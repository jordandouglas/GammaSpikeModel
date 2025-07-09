# Simulating trees and datasets using ReMASTER and the gamma spike model

## Dependencies

- BEAST v2.7, which could be directly called up upon by `beast`
- GammaSpikeModel package installed
- Artifact feastMCMC.jar built from this repository
- R v4, available from command line as `Rscript` 

## Instructions

1. Simulate 200 FBD trees and simulate 200 datasets of 500-bp long under these trees (with the gamma spike model),  

```bash prepare.sh --nsims 200 --seqlength 500```

2. Run parallel jobs with one cpu per replicate using SLURM,

```bash submit_all_runs.sh --nsims 200```

## XML files

- `truth.xml` - simulating the FBD tree and other parameters. The stubs are removed from the tree and the number of stubs per branch are reported.
- `seqsim.xml` - simulates sequences down the tree, which is specified in json format and parsed into the xml file through the `-df` option of beast.
- `run.xml` - performs MCMC on a dataset specified by the `-df` option. For mixing convenience, the initial state is the true state.


