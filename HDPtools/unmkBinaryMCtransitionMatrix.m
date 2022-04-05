function pvec = unmkBinaryMCtransitionMatrix(T)
% pvec = unmkBinaryMCtransitionMatrix(T)
%
% un-Builds sparse transition matrix for binary Markov chains, where number of
% states is given by 2^N, where N is the Markov order.
%
% Output: pvec [N x 1] - vector of probabilities of "1" given each
%        state, where states are assumed to be sorted in increasing number 
%       
% Input: T - [N x N] - sparse matrix, whose i'th row contains
%         the probability of transitioning to each state from state i.
%

N = size(T,1);  % number of states
Nh = N/2;  % half the number of states

% Compute indices for sparsity pattern
I = (1:N)';
J = reshape([Nh+1:N;Nh+1:N],[],1);
inds = sub2ind([N,N],I,J); % indices

% remove relevant states from sparse matrix
pvec = full(T(inds));