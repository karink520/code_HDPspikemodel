function v = cellcelleval(funptr, cellarray)
%  v = cellcelleval(funptr, cellarray)
%
%  evaluates a function 'fun' on each element of a cell array, 
%  returns a cell array of the same size
%
%  Inputs: funptr - must return a scalar
%          cellarray - cell array whose cells can be evaled by funptr
%
%  examples:  
%    cellmaxes = celleval('max', cellarray)
%    absmaxes  = celleval(@(x)max(abs(x(:))), cellarray);

csize = size(cellarray);
v = cell(csize);
for j = 1:csize(1)
    for k = 1:csize(2)
    v{j,k} = feval(funptr, cellarray{j,k});
    end
end
