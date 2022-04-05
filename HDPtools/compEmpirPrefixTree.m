function SS = compEmpirPrefixTree(spktrain,depth)
% SS = compEmpirPrefixTree(spktrain,depth)
%
% Compute Empirical Prefix Tree probabilities p(1|context) from spike train
%
% (Note: this is an inefficient way to compute it; Not something we should
% need in general)



SS = cell(depth+1,1);

SS{1} = mean(spktrain); 

for jj=2:depth+1
    X = extractMCstatesFromSpikes(spktrain,jj-1);
    [~,Tempir] = countMCstateTransitions(X,jj-1);
    pvec = unmkBinaryMCtransitionMatrix(Tempir); % estimated transition probabilities
    SS{jj} = reshape(pvec,2,[]);
end

