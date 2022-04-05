function pvecs = sampleBinaryHDPparams(pmu,alphas,nlevels)
% pvecs = sampleBinaryHDPparams(alphas,nlevels)
%
% Sample Beta RVs hierarchically, as drawn from a hierarchical DP, with
% concentration parameters 'alphas'
%
% Input:     pmu [1 x 1] - mean (avg probability of a spike in a bin)
%         alphas [nlevels+1 x 1] - concentration params (ok if scalar)
%
% Output: pvecs [nlevels+1 x 1] - cell array of p(1|context) for all levels
%         Each element is a [2 x level^2] matrix, with each row drawn from
%         identical distribution (formed by taking the pvec from the higher
%         level, reshaping it into a column vector and then taking
%         transpose).

% Parse inputs
if nargin > 2
    if length(alphas)==1
	alphas = repmat(alphas,nlevels+1,1);
    end
else
    nlevels = length(alphas)-1;
end
pvecs = cell(nlevels+1,1); % create cell array
 
% Draw top-level RV:  p(spike|empty context)
pvecs{1} = betarnd(alphas(1)*pmu,alphas(1)*(1-pmu));

% Draw parameters for each successive level
for j = 1:nlevels
    pprev = pvecs{j}(:)';
    pvecs{j+1} = betarnd(...
	alphas(j+1)*repmat(pprev,2,1), ...
	alphas(j+1)*repmat(1-pprev,2,1));
end
