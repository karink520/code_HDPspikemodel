function [ g ] = wordprob( data, K, alphas,pmu )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% K = length of context
% alphas size (K+1,1)
% here the w is 0

%  Output:  g, a 2^K x 1 vector of transition probabilities to 0 from 000,001,010,
%  (or whatever, but in that order)

if K == 0  %no context    (mix of data and prior)
    cu = length(data);
    cuw = (count(data));
    cuw = cuw(1);
    g = cuw/(alphas(1) + cu) + alphas(1)/(alphas(1) + cu) * pmu;

else 
    %c = cell{K + 1, n};

   


%cuw = # of u followed by a w
%cu = # appearances of u 
%theta is the concentration parameter
theta = alphas(K+1)*ones(2^K,1);

    cuw = countBlocks2(K, data); %counts all the blocks of length K+1
    indsEven = [2:2: length(cuw)];
    indsOdd = [1:2:length(cuw)];
    cu = cuw(indsEven) + cuw(indsOdd);  %extracts total of each length-K context
    cuw = cuw(indsOdd);   %extracts the one that end in 0
     
    g = cuw./(theta + cu) + (theta./(theta + cu)) .* repmat(wordprob(data, K-1, alphas, pmu),2,1);
end

end

