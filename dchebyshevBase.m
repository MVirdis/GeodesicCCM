function ys = dchebyshevBase(degs, xs)
%DCHEBYSHEVBASE(degs, xs) Returns the derivative of a chebyshev polynomial
%(first type) base.
%   Inputs: degs => an array of integers listing all the degrees of the
%   Chebyshev polynomials. For example, degs=0:N, where N is the maximum
%   degree.
%           xs {1,num_cpnts} => points at which to evaluate the polynomials
%   Outputs: ys {N+1,num_cpnts} => matrix of polynomial derivative
%   evaluations. Along the rows all the varying increasing polynomial
%   degrees. Along the columns the evaluation at each different point.

if isa(xs, 'double')
    ys = zeros(length(degs), length(xs));
elseif isa(xs, 'sym')
    ys = sym( zeros(length(degs), length(xs)) );
end

for i=2:length(degs)
    ys(i,:) = (i-1)*chebyshevU(i-2, xs);
end

end