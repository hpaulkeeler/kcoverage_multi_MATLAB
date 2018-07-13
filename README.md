# kcoverage_multi_MATLAB

If you use this code for published work, please cite paper[2] listed below.

The signal-to-interference-plus-noise ratio or SINR is considered a key performance metric in wireless communications based on information theoretic arguments. These Matlab scripts calculate (via integration or simulation) the SINR-based k-coverage probability of a typical user (or UE) in a multi-tier network with constant base station density (for each tier) and a common path-loss exponent based on Poisson network models outlined in [1] (single-tier) and [2] (multi-tier). The results hold for the whole SINR domain ie when the SINR threshold tau>0, and not just tau>1, and for arbitrary fading/shadowing or "propagation effects".

The calculations can take markedly longer for small tau values eg less than 0.1 (or -10 in dB).

Files that start with "Test" can be executed without passing parameters eg open file TestIntSTVsMT.m and hit F5. Files that start with "fun" are functions files that require parameters to be passed to them. 

The expression for the k-coverage coverage probability consists of a product of two integrals. The first integral (over the real line), calculated in the file funIn.m, captures the noise, density and propagation effects of the network model. For a so-called "interference limited network" where the noise is set to zero (ie W=0), this integral reduces to a closed-form solution (equation (15) in [2]). If noise is present in the network, the numerical integration method implemented in funIn.m breaks down for large values of n (~ 20, corresponding to very small values of SINR threshold eg tau <-11 dB)

The second integral, calculated respectively in funJn.m and funJnMT.m for single-tier and multi-tier networks, is defined over a (k-1)-dimensional hyper cube. The integral is calculated (or estimated) via a quasi-Monte Carlo scheme using Sobol points, which theoretically outperform regular Monte Carlo methods. Quasi-random points such as Sobol points can be chosen based on a thorough analysis of the integral kernel. However, this has not been done here as Sobol points, which perform generally well for "regular" integrals, appear to already give fast and accurate solutions.

These results hold for any fading and/or shadowing probability distribution as long as they adhere to a moment condition, which isdue to a property, sometimes called propagation invariance[3], based on using a Poisson process for the base station configuration, that says that propagation effects only depend one key moment, see [2] and [3] for more details on propagation invariance and the equivalence of Poisson networks.

The multi-tier work uses factorial moment measures of the SINR process and the closely related STINR process (knowing one is equivalent to knowing the other via simple algebraic relation, Corollary 8 in [2]). Of course, setting the number of tiers to one in mult-tier files gives single-tier results.

In the case of no noise (ie W=0), single-tier networks are scale-invariant (ie independent of base station density, lambda). For multi-tier networks the ratios between base station densities (and not their absolute values) matter when there is no noise. See [4] for more details on the multi-tier network model (but only tau>1 case is covered).

File descriptions:

TestIntSTProbCovk.m: Calculates k-coverage probability for a single-tier cellular network with log-normal shadowing based on a Poisson model outlined in [1].

TestIntSTVsMT.m: Calculates the k-coverage probability in a multi-tier network and a single-tier equivalent network. The single-tier equivalent has the same propagation processes with random SINR threshold values averaged; see [3] for more details on equivalent Poisson networks.

TestSimVsIntMT.m: Calculates the k-coverage probability in a multi-tier network by both simulating a network model and using an integration method outlined [2].

TestDeltaProbIntExamples.m: Calculates (via integration) the Delta probability for 2-base station signal combination and second-strongest signal removed in a multi-tier cellular network based on a Poisson model outlined in [2].

TestDeltaProbSimExamples.m: Estimates (via simulation) the Delta probability in the setting of TestDeltaProbIntExamples.m.

Author: H.P. Keeler, Inria Paris/ENS, 2014

References: 
[1] H.P. Keeler, B. Błaszczyszyn and M. Karray, 
'SINR-based coverage probability in cellular networks with arbitrary shadowing', ISIT, 2013 
[2] B. Błaszczyszyn and H.P. Keeler, 'Studying the SINR process in Poisson networks by using its factorial moment measures', 2015, published in IEEE Transactions on Information Theory. 
[3] Błaszczyszyn and H.P. Keeler, 'Equivalence and comparison of random heterogeneous networks', presented at WDN-CN, 2013 
[4] Dhillon, Ganti, Adnrews, Baccelli, 'Modeling and analysis of K-tier downlink heterogeneous cellular networks', JSAC, 2012
