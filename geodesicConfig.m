% Filename: geodesicConfig.m
% Author: Mario Virdis
% Created: 2025-04-16
% Description: See function header.
% Version: 1.0

function config = geodesicConfig(N,num_cpnts,nx,W,dW)
%geodesicConfig Returns the optimization parameters given a configuration
%structure config, with fields N (max degree of Chebyshev polynomials) and
%num_cpnts (number of collocation points), W the inverse contraction
%metric, and dW the derivative of the inverse CM, nx the dimensions of the
%state. Both W(x) and dW(x) must be functions. The returned structure
%contains all the ingredients to solve the geodesic optimization problem.

% Init CCQ scheme
[cpnts, weights] = ingredientsCCQ(0,1,num_cpnts);

% Compute basis polynomials
basis = chebyshevBase(0:N, cpnts);   % {N+1, num_cpnts}
dbasis = dchebyshevBase(0:N, cpnts); % {N+1, num_cpnts}

config.N = N;
config.num_cpnts = num_cpnts;
config.nx = nx;
config.W = W;
config.dW = dW;
config.cpnts = cpnts;
config.weights = weights;
config.basis = basis;
config.dbasis = dbasis;

end
