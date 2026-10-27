# Exact cumulant identifiability for gamma-spike burst sums

For B~Poisson(beta), independent G_i~Gamma(shape alpha,scale theta), and D=sum(i=1,...,B+1)G_i, this packet studies the finite cumulant map on all alpha,theta,beta>0. Full-distribution identifiability is already known; the contribution here is a sharp finite-information classification.

- Every positive three-cumulant fiber has an explicit admissible-quartic inverse and contains at most two parameter triples.
- On variance/mean²<=1/6 the three-cumulant map is injective. For each fixed ratio above 1/6, some fibers contain two triples; this does not classify every such point as ambiguous.
- Four cumulants globally identify the positive triple and are the least global order. The proof excludes all possible three-cumulant collisions analytically.

For example, (alpha,theta,beta)=(90/67,67/360,3) and (240/109,109/720,2) share cumulants(1,269/720,57241/259200), while their fourth cumulants differ by40043/7464960. These are exact rational parameters, not fitted data.

Install the single pinned dependency and replay from the repository root:

```sh
python -m pip install -r examples/exact-cumulants/requirements.txt
python examples/exact-cumulants/reproduce.py
python -O examples/exact-cumulants/reproduce.py
```

For the exact shared three-cumulant vector, the helper returns both candidates; including the first witness's fourth cumulant selects that triple:

```sh
python examples/exact-cumulants/gamma_cumulants.py 1 269/720 57241/259200 --json
python examples/exact-cumulants/gamma_cumulants.py 1 269/720 57241/259200 5417707/31104000 --json
```

The helper accepts three or four exact rational cumulants. It uses rational polynomial gcds and Sturm root isolation to return all positive candidates, preserving exact algebraic roots. Floating inputs and noisy moment fitting are outside its contract. The symbolic checker verifies displayed identities and the collision; the complete-domain root-count and sign conclusions require [proof.md](proof.md).

In the GammaSpikeModel notation, alpha=spikeShape and alpha*theta=spikeMean. This theorem concerns the homogeneous burst component with the stated Poisson count law. Full branch distance includes rate*length+spike, sampled-ancestor branches can follow other conventions, and pooled heterogeneous branches need not have this distribution. The result is not a complete BEAST identifiability theorem, a claim of unstable or incorrect MCMC, or an empirical estimator.

[sources.md](sources.md) credits Dieselhorst–Stadler's known full-law result and Douglas et al.'s original model/software, and records the exact B+1 convention. Development and independent checking use substantial AI assistance. Review status appears in status.json; human verification and priority are not established. This is an original code/proof packet, with no journal manuscript or submission.

The corresponding independently replayed lab project is [AI Biology Lab PR37](https://github.com/ArtificialZeng/AI-Biology-Lab/pull/37). This optional example uses SymPy separately from the BEAST Java package.
