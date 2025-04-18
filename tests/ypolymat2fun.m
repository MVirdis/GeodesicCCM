function ff = ypolymat2fun(ypolymat_txt, ind_names)
%YPOLYMAT2FUN(ypolymat, ind_names) Convert Yalmip Polynomial Matrix to
%MATLAB function
%   Inputs: ypolymat_txt => the Yalmip polymat obtained by calling
%sdisplay on the optimization variable. Typically a cell array
%           ind_names => names of the indipendent variables as a cell array
%           of strings
%   Output: ff => an anonymous function in one input

polystrs = ypolymat_txt;
for i=1:length(ind_names)
    polystrs = replace(polystrs, ind_names{i}, ['x(' num2str(i) ')']);
end
funcode = '@(x) ([';
nrows = size(ypolymat_txt,1);
ncols = size(ypolymat_txt,2);
for irow=1:nrows
    for icol=1:ncols
        if icol < ncols
            funcode = strcat(funcode, [polystrs{irow,icol} ',']);
        else
            funcode = strcat(funcode, polystrs{irow,icol});
        end
    end
    if irow < nrows
        funcode = strcat(funcode, ';');
    else
        funcode = strcat(funcode, '])');
    end
end

ff = str2func(funcode);

end
