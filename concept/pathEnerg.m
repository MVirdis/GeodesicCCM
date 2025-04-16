function [E, dE, Es_, dEs_] = pathEnerg(W,dW,gamma,dgamma,basis,dbasis,w)
%PATHENERG computes the approximated energy along a path gamma
%evaluated at a set number of points, weighted by weights
%w_til, using the (inverse) Riemannian metric W_til evaluated along the
%path at the same locations. Returns a scalar. W is a function.

num_cpnts = size(gamma,2);
nx = size(gamma,1);
N = size(basis,1)-1;
if length(w) ~= num_cpnts
    error('The number of weights must be the same as the number of points at which the gradients dgamma/ds are evaluated');
end

Es = zeros(1,num_cpnts);
dEs = zeros(nx,N+1,num_cpnts);
for k=1:num_cpnts
    gamma_k = gamma(:,k); % {nx, 1}
    dgamma_k = dgamma(:,k); % {nx, 1}
    W_k = W(gamma_k); % {nx,nx} - the value of the inverse metric at the current state
    M_k = W_k\eye(nx); % {nx,nx} - fast inverse
    dW_k = dW(gamma_k); % {nx,nx} - the value of the inverse metric derivative at the current state
    Es(1,k) = dgamma_k'*M_k*dgamma_k*w(k);
    for i=1:nx
        for j=1:N+1
            % Compute dE/dc_ij
            dEs_1 = 2*dgamma_k'*M_k(:,i)*dbasis(j); % {1,1}
            dEs_2 = -dgamma_k'*M_k*dW_k*M_k*dgamma_k*basis(j); % {1,1}
            dEs(i,j,k) = (dEs_1 + dEs_2)*w(k);
        end
    end
end

E = sum(Es,2);
dE = sum(dEs,3);

if nargout > 2
    Es_ = Es;
    dEs_ = dEs;
end

end