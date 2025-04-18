% Filename: evalPath.m
% Author: Mario Virdis
% Created: 2025-04-16
% Description: See function header.
% Version: 1.0

function [xs, dxs] = evalPath(config, c, pnts, recompute)
%PATH evaluates a parametrized curve on a sequence of points of a virtual
%coordinate s. Can return (in addition) the speed along the curve.
%   Detailed explanation goes here

N = config.N;
nx = config.nx;

if nargin == 3
    recompute = false;
end

% Compute basis polynomials
if recompute
    basis = chebyshevBase(0:N, pnts);   % {N+1, num_cpnts}
else
    basis = config.basis;
end

% Decision matrix from the vector
C = reshape(c, [N+1, nx])';
xs = C*basis;

if nargout > 1
    if recompute
        dbasis = dchebyshevBase(0:N, pnts); % {N+1, num_cpnts}
    else
        dbasis = config.dbasis;
    end
    dxs = C*dbasis;
end

end
