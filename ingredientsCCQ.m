% Filename: ingredientsCCQ.m
% Author: Mario Virdis
% Created: 2025-04-16
% Description: See function header.
% Version: 1.0

function [cpnts, weights] = ingredientsCCQ(a,b,num)
%INGREDIENTSCCQ(a,b,num) returns collocation points and weights for the
%Clenshaw–Curtis quadrature scheme.
%   Inputs: a, b => the interval [a,b] on which the integral is defined. (a
%must be less than b)
%           num => number of collocation points.
%   Outputs: cpnts => collocation points in the interval [-1,1] {1,num}
%            weights => the weights for the integrand {1,num}
%   We assume that, because we are using chebyshev polynomials, the points
%   have to lie in [-1,1].

if a >= b
    error('a must be strictly less than b');
end

k1 = 0.5*(b-a);
% k2 = 0.5*(a+b);

cpnts = cos(pi*(0:num-1)/(num-1)); % {1,num}
weights = weightsCCQ(num-1);       % {1,num}

% Scale to interval [a,b]
weights = k1*weights;

end

function w = weightsCCQ(N)
%Code edited from the accompanying code to arXiv:1607.04340v2 by K. Leung
    w = zeros(1,N+1);
    if even(N)
        w(1) = 1/(N^2 - 1);
        w(end) = w(1);
        for s=1:floor(N/2)
            wi = 0;
            for j=0:N/2
                if j == 0
                    wi = wi + 0.5 /(1-4*j^2) * cos(2*pi*j*s/N);
                elseif j == N/2
                    wi = wi + 0.5 /(1-4*j^2) * cos(2*pi*j*s/N);
                else 
                    wi = wi + 1 /(1-4*j^2) * cos(2*pi*j*s/N);
                end
            end
            w(s+1) = 4/N*wi;
            w(N-s+1) = w(s+1);
        end
    else
        w(1) = 1/N^2;
        w(end) = w(1);
        for s = 1:floor((N-1)/2)
            wi = 0;
            for j = 0:(N-1)/2
                if j == 0
                    wi = wi + 0.5 /(1-4*j^2) * cos(2*pi*j*s/N);
                elseif j == N/2
                    wi = wi + 0.5 /(1-4*j^2) * cos(2*pi*j*s/N);
                else 
                    wi = wi + 1 /(1-4*j^2) * cos(2*pi*j*s/N);
                end
            end
            w(s+1) = 4/N*wi;
            w(N-s+1) = 4/N*wi;
        end
    end
end
