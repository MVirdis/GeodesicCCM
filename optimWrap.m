% Filename: optimWrap.m
% Author: Mario Virdis
% Created: 2025-04-16
% Description: A wrapper to pathEnergy for use with optimizers.
% Version: 1.0

function [E,dE] = optimWrap(config, c, grad)
%TODO doc

if grad
    dataout = pathEnergy(config, c);
    E = dataout.E;
    dE = dataout.dE;
else
    dataout = pathEnergy(config, c, {'E'});
    E = dataout.E;
end

end
