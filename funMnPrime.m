% funMnPrime calculates and returns the n-th moment measure of the STINR
% process defined by equation (8) and (9) in [2]
% Mn=funMn(tValues,betaConst,n,x,gammaConst,numbMC)
% Mn is the moment measure of the STINR process Z'
% tValues is the n-valued of STINR threshold values
% betaConstant is path-loss exponent
% n is integer parameter
% betaConst, n, and numbMC are scalars. T_n is a vector of length n
% numbMc is number of sample points for  quasi-MC integration
%
% Author: H.P. Keeler, Inria/ENS, Paris and University of Melbourne, 2018
%
% References
% [1] H.P. Keeler, B. Błaszczyszyn and M. Karray,
% 'SINR-based coverage probability in cellular networks with arbitrary
% shadowing', ISIT, 2013
% [2] B. Błaszczyszyn  and H.P. Keeler, 'Studying the SINR process in
% Poisson networks by using its factorial moment measures', IEEE TOIT,
% 2015
%
% NOTE: If you use this code in published work, please cite paper[2] by
% Błaszczyszyn and Keeler, as listed above.


function Mn=funMnPrime(tPrimeValues,betaConst,n,x,gammaConst,numbMC)

tPrimeValues=tPrimeValues(:); %converts tValues to a column vector


if (gammaConst*sum(tPrimeValues))<1; %eq. (22) in [2]
    In=funIn(betaConst,n,x); %eq. (12) in [1] or eq. (13) in [2]
    t_hat=gammaConst*tPrimeValues./(1-gammaConst*sum(tPrimeValues));
    JnMT=funJnMT(t_hat,betaConst,n,numbMC); %eq. (16) in [2]
    Mn=prod(t_hat.^(-2/betaConst)).*In.*JnMT.*factorial(n); %eq. (21) in [2]
else
    Mn=0;
end
