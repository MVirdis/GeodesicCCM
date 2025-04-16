function [E,dE] = energy_fun(W,dW,basis,dbasis,w,nx, cx)
%ENERGY_FUN Returns path energy given decision variables. To be used with
%fmincon or other similar optimizers. Rebuilds the basis given the
%coefficients and returns the path energy and its gradient.

N = size(basis,1)-1;

% Rebuild coefficient matrix
C = reshape(cx, [N+1,nx])';

gamma = C*basis;
dgamma = C*dbasis;

[E,dE] = pathEnerg(W,dW,gamma,dgamma,basis,dbasis,w);

end
