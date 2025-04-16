% Filename: pathEnergy.m
% Author: Mario Virdis
% Created: 2025-04-16
% Description: See function header.
% Version: 1.0

function dataout = pathEnergy(config, c, request)
%PATHENERGY computes the approximated energy along a path gamma
%evaluated at a set number of points, weighted by weights
%w_til, using the (inverse) Riemannian metric W_til evaluated along the
%path at the same locations. Returns a scalar. W is a function.
%It can optionally accept a third input a cell array with strings:
%{'E','dE','Es','dEs'} (or a subset of those). Made for efficiency.

% Unpack parameters
W = config.W;
dW = config.dW;
basis = config.basis;
dbasis = config.dbasis;
w = config.weights;
num_cpnts = config.num_cpnts;
nx = config.nx;
N = config.N;

dataout = struct();

% Init flags
Ef = false;
dEf = false;
Esf = false;
dEsf = false;

if nargin == 2
    % Default to energy and gradient, without arrays
    Ef = true;
    dEf = true;
    Esf = false;
    dEsf = false;
else
    for i=1:length(request)
        switch request{i}
            case 'E'
                Ef = true;
            case 'dE'
                dEf = true;
            case 'Es'
                Esf = true;
            case 'dEs'
                dEsf = true;
        end
    end
end

% Decision matrix from the vector
C = reshape(c, [N+1, nx])';
gamma = C*basis;
dgamma = C*dbasis;

Es = zeros(1,num_cpnts);
dEs = zeros(nx,N+1,num_cpnts);
if dEf
    dEsw = zeros(nx,N+1,num_cpnts);
end
for k=1:num_cpnts
    gamma_k = gamma(:,k); % {nx, 1}
    dgamma_k = dgamma(:,k); % {nx, 1}
    W_k = W(gamma_k); % {nx,nx} - the value of the inverse metric at the current state
    M_k = W_k\eye(nx); % {nx,nx} - fast inverse
    Es(1,k) = dgamma_k'*M_k*dgamma_k;

    if dEf || dEsf % - avoids to compute gradient info if not requested
        dW_k = dW(gamma_k); % {nx,nx} - the value of the inverse metric derivative at the current state
        for i=1:nx
            for j=1:N+1
                % Compute dE/dc_ij
                dEs_1 = 2*dgamma_k'*M_k(:,i)*dbasis(j); % {1,1}
                dEs_2 = -dgamma_k'*M_k*dW_k*M_k*dgamma_k*basis(j); % {1,1}
                dEs(i,j,k) = dEs_1 + dEs_2;
                dEsw(i,j,k) = (dEs_1 + dEs_2)*w(k);
            end
        end
    end
end

if Ef
    E = sum(Es.*w,2);
    dataout.E = E;
end
if dEf
    dE = sum(dEsw,3);
    dataout.dE = dE;
end
if Esf
    dataout.Es = Es;
end
if dEsf
    dataout.dEs = dEs;
end

end
