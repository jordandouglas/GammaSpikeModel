
# Simulating trees and datasets under the gamma spike model


## Dependencies

- BEAST 2.7 - assumed to be installed in the `~/beast/bin/` directory
- GammaSpikeModel package installed
- R v4, available from command line as `Rscript` 


## Instructions



1. Simulate 100 FBD trees and then simulate 100 datasets under these trees (with gamma spike model),  


```bash prepare.sh```





2. Then to run MCMC on a single replicate:

```
cd templates
cd rep1
~/beast/bin/beast -df var.seq.json ../../run.xml
```

3. Or, to run MCMC on all replicates one at a time:

```
bash run.sh
```



4. Or, you can run MCMC on all replicates in parallel across multiple threads (default 4). **Warning**: these jobs will run in the background, so do not run this script unless you are prepared to monitor/kill the jobs


```
bash runParallel.sh 
```


## XML files

- `truth.xml` - simulating the FBD tree and other parameters. The stubs are removed from the tree and the number of stubs per branch are reported. By default, there are 40 extant taxa and a variable number of non-extant taxa.

- `seqsim.xml` - simulates sequences down the tree, which is specified in json format and parsed into the xml file through the `-df` option of beast.

- `run.xml` - performs MCMC on a dataset specified by the `-df` option. For mixing convenience, the initial state is the true state.


