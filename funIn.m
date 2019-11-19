% funIn calcualtes and returns the In integral (eq (12) in [1] or eq (13)
% in [2])
% In=funIn(betaConst,n,x)
% In = I_{n,\beta} integral (scalar) value
% betaConstant = path-loss exponent
% n = integer parameter
% x = variable that incorporates model parameters
% That is, x=W*a^(-2/betaConst) where a is given by eq (6) in [1]
% betaConst, n and x are scalars
%
% WARNING: for x>0 and high values of n (corresponding to low values of SINR
% threshold eg -11 dB) the integration method can fail (returning NaN). This is due
% to a singularity in the integral kernel. A couple different change of variables
% have been tried with different (Matlab) quadratuture methods. The best
% performing combination has been chosen. No doubt better numerical methods
% exist (possibly combining an asymptotic expansion at the singularity).
%
% Author: H.P. Keeler, Inria Paris/ENS, 2014
%
% References:
% [1] H.P. Keeler, B. Blaszczyszyn and M. Karray,
% 'SINR-based coverage probability in cellular networks with arbitrary
% shadowing', ISIT, 2013
%
% NOTE: If you use this code in published work, please cite paper[1] by
% Keeler, Blaszczyszyn and Karray, as listed above.



function In=funIn(betaConst,n,x)
% Calculates I_n with numerical integration or analytic solution (if x=0 
% (ie W=0))x <>0 uses integral (which can handle singularities). An older 
% version of the code used the integration function quadgk.
% https://mathworks.com/help/matlab/ref/integral.html
% x=0 uses analytic solution

C=gamma(1+2./betaConst).*gamma(1-2./betaConst); %constant C';eq. (13) in [1] 
if x==0
    %analytic solution
    In=(2./betaConst).^(n-1)./C.^n; %eq. (14) in [1]
else
    %numerical solution         
    %latest quadrature of (12) in [1] after change of variable
    %u=t./(1-t); change of variable
    ft=@(t)t.^(2*n-1).*(1./(1-t)).^(2*n+1).*exp(-(t./(1-t)).^2-(t./(1-t)).^(betaConst).*x.*((gamma(1-2./betaConst)).^(-betaConst/2)));
    In=2^n./betaConst.^(n-1).*integral(ft,0,1)./C.^n/factorial(n-1);      
    
end

