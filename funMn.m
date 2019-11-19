% funMn calculates and returns the n-th moment measure of the SINR
% process defined by equation (5) and (6) in [2]
% Mn=funMn(tValues,betaConst,n,x,gammaConst,numbMC)
% Mn is the moment measure of the SINR process Z
% tValues is the n-valued of SINR threshold values
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

function Mn=funMn(tValues,betaConst,n,x,gammaConst,numbMC)
tValues=tValues(:); %converts tValues to a column vector
tPrimeValues=tValues./(1+gammaConst*tValues);

if (gammaConst*sum(tPrimeValues))<1; %eq. (22) in [2]
    In=funIn(betaConst,n,x); %eq. (12) in [1] or eq. (13) in [2]
    t_hat_prime=gammaConst*tPrimeValues./(1-gammaConst*sum(tPrimeValues));
    t_hat=t_hat_prime./(1+gammaConst*t_hat_prime); %calculate t_hat prime values
    JnMT=funJnMT(t_hat,betaConst,n,numbMC); %eq. (16) in [2]
    Mn=prod(t_hat.^(-2/betaConst)).*In.*JnMT.*factorial(n); %eq. (21) in [2]

else
    Mn=0;
end
