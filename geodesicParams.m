% Filename: geodesicParams.m
% Author: Mario Virdis
% Created: 2025-04-16
% Description: See function header.
% Version: 1.0

function params = geodesicParams(config, xeq, x)
%GEODESICPARAMS(config, c, xeq, x) Summary of this function goes here
%   Detailed explanation goes here

nx = config.nx;
N = config.N;

C0 = [0.5*(xeq+x), 0.5*(x-xeq), zeros(nx,N-1)]; % Initial guess: straight line connecting xeq and x
% Flatten (! careful MATLAB uses column-major)
C0_T = C0';
c0 = C0_T(:);

% Boundary conditions
%    ! this formulation assumes that the decision matrix C is flattened to
%    a decision vector in row-major order. (MATLAB flattens column-major,
%    so must transpose first)
base1 = chebyshevBase(0:N, -1)';
base2 = chebyshevBase(0:N, 1)';
Aeq1 = base1;
Aeq2 = base2;
for i=1:nx-1
    Aeq1 = blkdiag(Aeq1, base1);
    Aeq2 = blkdiag(Aeq2, base2);
end
Aeq = [Aeq1; Aeq2];
beq = [xeq;x];

params.C0 = C0;
params.c0 = c0;
params.Aeq = Aeq;
params.beq = beq;

% Cost function
params.f = @(c) ( optimWrap(config, c, true) );

end
