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
