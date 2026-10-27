# Four cumulants identify the shifted compound-Poisson gamma model

## Model and scope

Let \(B\sim\operatorname{Poisson}(\beta)\) and let \(G_i\) be independent
\(\operatorname{Gamma}(\alpha,\theta)\) variables, independent of \(B\), with
shape \(\alpha>0\), scale \(\theta>0\), and \(\beta>0\). Set
\(D=\sum_{i=1}^{B+1}G_i\). The source's conditional gamma law and cumulant
generating function specify this \(B+1\) model. An inconsistent summation
index in its identifiability theorem is not used as a counterexample.

Write \(F_m(\alpha,\theta,\beta)=(\kappa_1,\ldots,\kappa_m)\).
The source already proves identifiability from the complete distribution.
Here the question concerns the finite cumulant map on the full positive
parameter domain. No empirical estimation or noise-stability assertion is made.

## Statements

1. Every positive fiber of \(F_3\) is obtained by the explicit admissible
   quartic inversion below and contains at most two triples. The bound is sharp.
2. Put \(c=\kappa_2/\kappa_1^2\). The restriction to \(c\le1/6\) is injective.
   For every fixed \(c>1/6\) and every positive prescribed mean, some
   three-cumulant fiber has two triples. This does not say that every
   parameter point with \(c>1/6\) is ambiguous.
3. \(F_4\) is injective on all \((0,\infty)^3\). Four is the least global
   cumulant order with this property.

### Complete inversion and root-count rule

For proposed real cumulants \(k_1,k_2,k_3\), if \(k_1\le0\) or \(k_2\le0\)
there is no positive fiber. Otherwise put
\[
 \mu=k_1,\qquad c=k_2/\mu^2,\qquad t=k_3/\mu^3,
 \qquad d_c(z)=z^2-z+c.
\]
Let
\[
 \mathcal D_c=\{0<z<1:d_c(z)>0\},\qquad
 f_c(z)=2c^2+cz^2-cz-z^4+z^3.
\]
For every real root \(z\in\mathcal D_c\) of
\[
 z^4-z^3-cz^2+cz+t-2c^2=0,                         \tag{1}
\]
there is exactly one triple in the fiber:
\[
 \theta=\mu d_c(z),\qquad
 \alpha=z/d_c(z),\qquad
 \beta=(1-z)/z.                                  \tag{2}
\]
There are no other triples.

Here is an exact count rule, including boundary and repeated-root cases.
The domain is \((0,1)\) if \(c>1/4\); for \(0<c\le1/4\) it is
\[
 (0,l_c)\ \cup\ (r_c,1),\qquad
 l_c=(1-\sqrt{1-4c})/2,\quad r_c=(1+\sqrt{1-4c})/2.
\]
At \(c=1/4\), the single point \(1/2\) is excluded.
There is one derivative root \(u_c\in(0,1/2)\); if \(c<1\) there is also
one \(v_c\in(1/2,1)\). They are the indicated roots of
\[
 g_c(z)=-4z^3+3z^2+2cz-c.                         \tag{3}
\]
Remove these derivative roots from \(\mathcal D_c\). On each remaining
connected open interval, count one solution exactly when \(t\) is strictly
between the two one-sided endpoint values of \(f_c\). Add one for each
removed derivative root that lies in \(\mathcal D_c\) and has \(f_c=t\).
Endpoints outside \(\mathcal D_c\) are never added. This rule and (2)
specify every positive fiber without a numerical parameter grid.

## Cumulants and the bijective normalization

The cumulant generating function is
\[
 K(s)=-\alpha\log(1-\theta s)
       +\beta((1-\theta s)^{-\alpha}-1).
\]
Consequently
\[
 \kappa_n=\theta^n\{\alpha(n-1)!+\beta(\alpha)_n\},\quad n\ge1.
\]
For \(h=\alpha\theta\), direct substitution gives
\[
 \kappa_{n+1}=(h+n\theta)\kappa_n
       -h^2(n-1)!\theta^{n-1}.                    \tag{4}
\]
In particular \(\mu=h(1+\beta)\). Thus \(z=h/\mu=1/(1+\beta)\)
lies in \((0,1)\). The \(n=1\) recurrence gives
\(\theta/\mu=d_c(z)>0\), and its converse is exactly (2).
Using \(n=2,3\) yields
\[
 \kappa_3/\mu^3=f_c(z),\qquad
 \kappa_4/\mu^4=j_c(z)
   =(z+3d_c(z))f_c(z)-2z^2d_c(z)^2.               \tag{5}
\]
This proves the inversion statement; all inequalities and denominators in
(2) have been retained.

## Fiber bound, sharp injectivity threshold, and completeness

On \((0,1/2)\), \(g'_c(z)=6z(1-2z)+2c>0\);
\(g_c(0)=-c<0\) and \(g_c(1/2)=1/4>0\). Hence there is exactly one \(u_c\).
The derivative \(g'_c\) has just one positive root,
\((3+\sqrt{9+24c})/12>1/2\). It follows that after \(1/2\), \(g_c\)
first increases and then decreases, when this turning point lies before one.
Since \(g_c(1)=c-1\), there is exactly one \(v_c\) if \(c<1\),
and none in \((1/2,1)\) if \(c\ge1\).
Thus \(f_c\) decreases, increases, and (only for \(c<1\)) decreases on the
successive intervals defined by these roots. This also proves the stated
monotone root-count rule, including critical values.

For \(c\ge1\), there is only the minimum \(u_c\), so every level has at
most two roots on \((0,1)\). For \(c<1\), use
\[
 f_c(z)-2c^2=z(1-z)(z^2-c),\qquad
 g_c(\sqrt c)=2c(1-\sqrt c)>0.
\]
Therefore \(u_c<\sqrt c<v_c\). A level below \(2c^2\) lies entirely left
of \(\sqrt c\), where there is just the one minimum; a level above it
lies entirely right of \(\sqrt c\), where there is just the one maximum.
The level \(2c^2\) has only \(z=\sqrt c\) in \((0,1)\).
Every level again has at most two roots. Restricting to \(\mathcal D_c\)
cannot increase that bound.

If \(0<c\le1/6\), evaluate at the two domain boundaries:
\[
 g_c(l_c)=l_c(6c-1)\le0,\qquad
 g_c(r_c)=r_c(6c-1)\le0.
\]
The derivative-root locations imply \(l_c\le u_c\) and \(r_c\ge v_c\).
Thus \(f_c\) is strictly decreasing on each admissible interval.
Also \(l_c<\sqrt c<r_c\), so their images lie strictly below and above
\(2c^2\), respectively, and cannot overlap. This proves injectivity,
including \(c=1/6\), where the derivative zeros are excluded boundaries.

For \(1/6<c<1/4\), the same evaluations are strictly positive, giving
\(u_c<l_c\) and \(r_c<v_c\). The minimum \(u_c\) is now an interior
point of \(\mathcal D_c\). For \(c\ge1/4\), it remains interior
(\(u_c<1/2\)). In each case a sufficiently small level above its strict
minimum has two distinct admissible roots. Scaling \(\mu>0\) in (2)
establishes the ambiguity assertion for every fixed \(c>1/6\).

## Four cumulants rule out all three-cumulant collisions

Suppose distinct \(z,w\in(0,1)\) have the same \(c,f_c\). Put
\(S=z+w\), \(T=zw\). Dividing \(f_c(z)-f_c(w)\) by \(z-w\) gives
\[
 (S-1)(c-S^2)+(2S-1)T=0.                         \tag{6}
\]
If \(S=1\), this is \(T=0\), impossible. Hence
\(c=S^2-(2S-1)T/(S-1)\). Substitution into the divided difference of
\(j_c\) gives the exact identity
\[
 \frac{j_c(z)-j_c(w)}{z-w}=-\frac{E(S,T)}{S-1},
\]
where
\[
 E=2S^3(S-1)^3-S(S-1)(8S^2-15S+6)T
        +(3S-2)(S-2)T^2.                        \tag{7}
\]
We show \(E<0\) on the entire relevant domain, even without using
the additional scale-positivity restrictions.

Always \(0<S<2\) and \(0<T<S^2/4\), because \(z,w\) are distinct.
For \(S>1\), additionally \(T>S-1\), since \((1-z)(1-w)>0\).
The endpoint identities are
\[
 E(S,0)=2S^3(S-1)^3,\quad
 E(S,S^2/4)=\frac{S^3(S-2)(3(S-1)^2+1)}{16}.
\]
For \(0<S\le2/3\), \(E\) is convex or linear in \(T\), and both
endpoint values are negative. Hence its entire intervening segment is negative.
For \(2/3<S<1\), its quadratic coefficient is negative and
\[
 \partial_T E(S,0)=-S(S-1)(8S^2-15S+6)<0.
\]
Indeed the convex quadratic \(8S^2-15S+6\) is negative throughout
\([2/3,1]\), as both endpoint values are negative. Thus \(E\) decreases
from the negative value at zero.

For \(1<S<2\), the quadratic coefficient is again negative and
\[
 E(S,S-1)=2(S-2)(S-1)^5<0,
\]
\[
 \partial_T E(S,S-1)=-(S-1)(8S^3-21S^2+22S-8)<0.
\]
The last cubic, with \(q=S-1>0\), is \(8q^3+3q^2+4q+1>0\).
Thus \(E\) decreases from a negative value throughout \(T>S-1\).
This exhausts all cases and proves that the fourth cumulants differ.
Together with (2), \(F_4\) is injective.

## An exact collision and minimality

In the order \((\alpha,\theta,\beta)\), take
\[
 A=(90/67,\ 67/360,\ 3),\qquad
 B=(240/109,\ 109/720,\ 2).
\]
Both have
\[
 (\kappa_1,\kappa_2,\kappa_3)
   =(1,\ 269/720,\ 57241/259200).
\]
Their fourth cumulants are respectively
\(5417707/31104000\) and \(31505167/186624000\), which differ.
All three parameters of both triples are strictly positive; even both
shape parameters exceed one. This proves noninjectivity for every order
\(m\le3\). The global \(F_4\) proof establishes minimality.

## Status and limitations

The proof is analytic and exact over the full positive parameter domain.
A finite symbolic checker may verify the displayed algebra and rational
witness, but does not by itself certify the root-count and sign arguments.
The result does not contradict the known full-law identifiability theorem,
establish stable inference from noisy empirical moments, or identify
biological histories from real data.

This is a proposed original code/proof contribution. Primary novelty/value
screening and fresh mathematical review are required; substantial AI
assistance is used. No human verification, author endorsement, priority
certification, or journal manuscript release is asserted.
