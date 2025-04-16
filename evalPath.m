function [xs, dxs] = evalPath(config, c, pnts)
%PATH evaluates a parametrized curve on a sequence of points of a virtual
%coordinate s. Can return (in addition) the speed along the curve.
%   Detailed explanation goes here

N = config.N;
nx = config.nx;

% Compute basis polynomials
basis = chebyshevBase(0:N, pnts);   % {N+1, num_cpnts}

% Decision matrix from the vector
C = reshape(c, [N+1, nx])';
xs = C*basis;

if nargout > 1
    dbasis = dchebyshevBase(0:N, pnts); % {N+1, num_cpnts}
    dxs = C*dbasis;
end

end
