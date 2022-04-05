function X = extractMCstatesFromSpikes(spks,depth)
% X = extractMCstatesFromSpikes(spks,depth)
%
% Extract Markov Chain states (of depth 'depth') from spike train.
%
% Inputs:
%    spks [N x 1] - spike train (vector of zeros and ones)
%   depth [1 x 1] - depth of Markov Chain
%
% Output:
% 

% Convert to column vector (if necessary)
if size(spks,2)>1
    spks = spks';
end

X = toeplitz([0;spks(1:end-1)],zeros(1,depth))*(2.^(depth-1:-1:0)')+1;

