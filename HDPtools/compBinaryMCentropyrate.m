function H = compBinaryMCentropyrate(p1trans,pstates)
% H = compBinaryMCentropyrate(p1trans,pstates)
%
% Computes entropy rate for binary markov chain from probabilities of state
% transitions and stationary distribution
%
% Inputs:
%   p1trans [N x 1] - for each state, probability of state ending in 1
%   pstates [N x 1] - asymptotic probability of each state
%
% Output:
%   H [1 x 1] - estimated entropy rate
%
% Formula:  H = \sum P(state) H(p(transitions|states))


ii = (p1trans>0)&(p1trans<1);

H = -sum(pstates(ii).*(p1trans(ii).*log2(p1trans(ii)) + ...
    (1-p1trans(ii)).*log2(1-p1trans(ii))));