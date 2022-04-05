function [T,Tpairs,ppairs] = mkBinaryMCtransitionMatrix(pvec)
% [T,Tpairs,ppairs] = mkBinaryMCtransitionMatrix(pvec)
%
% Builds sparse transition matrix for binary Markov chains, where number of
% states is given by 2^N, where N is the Markov order.
%
% Input: 
%     pvec [N x 1] - vector of p("1"| each state)
%                    (states ordered by HDP order)
%       
% Output: 
%        T [N x N] - sparse matrix, whose i'th row contains the probability
%                    of transitioning to each state from state i. 
%   Tpairs [N x 2] - list of 2 states that can be transitioned to for each
%                    of the N states
%   ppairs [N x 2] - probability of each of 2 states, for each state
%

N = length(pvec);  % number of states
Nh = N/2;  % half the number of states

% Compute indices for sparsity pattern
I = [(1:N)'; (1:N)'];
J = reshape([1:Nh;1:Nh],[],1); J = [J;J+Nh];

% Insert into sparse matrix
ppairs = [1-pvec,pvec];
T = sparse(I,J,ppairs(:));

if nargout > 1
    Tpairs = reshape(J,[],2);
end
    