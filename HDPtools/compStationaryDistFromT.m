function pmarg = compStationaryDistFromT(M)
% pmarg = compStationaryDistFromT(M)
%
% Computes left eigenvector of marginal state distribution (efficiently
% from a large, sparse matrix M).

SIZECUTOFF = 100;  % above this size, use "eigs" instead of "eig"

assert(size(M,1) == size(M,2), 'transition matrix must be square');

if size(M,1)<SIZECUTOFF
    % use FULL
    [pmarg,D] = eig(full(M'));
    [v,idx] = max(diag(D)); % eig doesn't always return sorted eigenvalues!!
    assert(v >= 0.99 && v <= 1.01, 'the maximum eigenvalue is not close to 1');
    pmarg = pmarg(:,idx)/sum(pmarg(:,idx)); % normalize to a probability vector
else
    % use SPARSE 
    [pmarg,v] = eigs(M', 1, 1.0001); % eigenvector with eigenvalue closest to 1
    assert(v >= 0.99 && v <= 1.01, 'the maximum eigenvalue is not close to 1!');
    pmarg = pmarg/sum(pmarg);  % make it normalize to 1
end
