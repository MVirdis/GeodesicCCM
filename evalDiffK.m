function [u,us_,pnts] = evalDiffK(config, c, dK, use_cpnts, Neval)
%EVALDIFFK evaluates a differential controller along a path.

if nargin == 3 || (nargin >= 4 && use_cpnts) % - reuse collocation points from the optimization process
    us = zeros(1,config.num_cpnts);
    [xs, dxs] = evalPath(config,c,config.cpnts);
    for i=1:config.num_cpnts
        xi = xs(:,i);
        us(1,i) = dK(xi)*dxs(:,i);
    end
    u = sum(us.*config.weights);

    if nargout > 1
        us_ = us;
    end
    if nargout > 2
        pnts = config.cpnts;
    end
elseif nargin == 5 % - re-approximate the integral with new collocation points
    us = zeros(1,Neval);
    [epnts, eweights] = ingredientsCCQ(0,1,Neval);
    [xs, dxs] = evalPath(config,c,epnts);
    for i=1:Neval
        xi = xs(:,i);
        us(1,i) = dK(xi)*dxs(:,i);
    end
    u = sum(us.*eweights);

    if nargout > 1
        us_ = us;
        pnts = epnts;
    end
else
    error('Wrong number of parameters');
end

end
