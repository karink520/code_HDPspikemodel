function [piSamps,topsmps] = gibbsSampleBinaryHDP(pmu,alphas,Tcounts,nsamps,npgrid,doplots)
% piSamps = gibbsSampleBinaryHDP(pmu,alphas,Tcounts,nsamps)
%
% Sample pi's at bottom level of an HDP model of binary spike train data.
%
% Input:  pmu [1 x 1] - mean at top level (avg probability of a spike in a bin)
%         alphas [N+1 x 1] - concentration params (can be scalar)
%         Tcounts [2^N x 2] - counts of transitions from each of 2^N states
%         nsamps [1 x 1] - number of Gibbs samples to return
%         npgrid [1 x 1] - number of bins to use for mesh on [0,1] for
%                          sampling probabilities pi
%         doplots [0 or 1] - boolean for whether or not to make plots
%
% Output: piSamps - Gibbs samples (transition probabilities)


nstates = size(Tcounts,1);
nlevels = log2(nstates);
piSamps = zeros(nstates,nsamps); % output variable for Gibbs samps
N0 = reshape(Tcounts(:,1),2,[]);  % counts of transitions to "0"
N1 = reshape(Tcounts(:,2),2,[]);  % counts of transitions to "1"

% Set bounds for p samples 
pmin = 1e-12;
pmax = 1-1e-12;

if length(alphas)==1
    alphas = repmat(alphas,nlevels+1,1);
end

% Make grid for sampling pi's  (option 1)
dp = 1/npgrid;
pgrid = (dp/2:dp:.99)';

% % Make grid for sampling pi's  (option 2)
% gmn = 1e-8;%min(pgrid);
% gmx = 1-gmn;%max(pgrid);
% dp = (gmx-gmn)/(npgrid-1);
% pgrid = (gmn:dp:gmx)';

% compute stuff
mgrid = 1-pgrid;
logpg = log(pgrid);
logmg = log(mgrid);

% Make top prior density
logpritop = log(betapdf(pgrid,alphas(1)*pmu,alphas(1)*(1-pmu)));

% Initialize with a sample from the prior
ppcell = sampleBinaryHDPparams(pmu,alphas,nlevels);

% ad hoc step: move each sample halfway back towards mean
ppcell = cellcelleval(@(x)((x+pmu)/2),ppcell); 

% now sample bottom level params using data (from Beta)
pup = repmat(ppcell{nlevels}(:)',2,1);
ppcell{nlevels+1} = max(pmin,min(pmax,betarnd(...
    alphas(nlevels+1)*pup+N1, ...
    alphas(nlevels+1)*(1-pup)+N0)));

% % Flag for making plots
% doplots = 0;

% Do Gibbs Sampling
for jj = 1:nsamps
    
    % sample top level param
    pdwn = ppcell{2}; % pis from below
    adwn = alphas(2); % alpha from below
    
    pdensity = logpritop + ...  % top beta prior 
	(adwn*pgrid-1).*(sum(log(pdwn))) + ...     % logli term 1
	(adwn*(mgrid)-1).*(sum(log(1-pdwn))) - ... % logli term 2
	2*betaln(adwn*pgrid,adwn*mgrid);  % logli normalizer
    pdensity = exp(pdensity-max(pdensity));
    
    
    ppcell{1} = randsample(pgrid,1,true,pdensity);
	
    if doplots
	subplot(211);
	plot(pgrid,pdensity);
	subplot(212);
	mkHDPtreeplot(ppcell);
	drawnow;
    end
    
    for jlevel = 1:nlevels-1;
	
	pup = reshape(repmat(ppcell{jlevel}(:)',2,1),1,[]); % pis from above
	aup = alphas(jlevel+1); % alpha from this level
	pdwn = ppcell{jlevel+2}; % pis from below
	adwn = alphas(jlevel+2); % alpha from below
    
	% Compute sampling density for pi's
	pdensity = ...
	    bsxfun(@times,logpg,(aup*pup-1))+ ... % beta pdf: 1st term
	    bsxfun(@times,logmg,(aup*(1-pup)-1))+ ... % beta pdf: 2nd term
	    bsxfun(@times,sum(log(pdwn),1),(adwn*pgrid-1))+ ... % beta li: 1st term
	    bsxfun(@minus,...
	    bsxfun(@times,sum(log(1-pdwn),1),(adwn*mgrid-1)), ...; % beta li: 2nd term
	    2*betaln(adwn*pgrid,adwn*mgrid)); % beta li: normalizer

	% Normalize density
	pdensity = exp(bsxfun(@minus,pdensity,max(pdensity,[],1)));
	
	% Sample each pi
	pcum = bsxfun(@rdivide,cumsum(pdensity),sum(pdensity)); % cumulative density
	rnds = rand(1,size(pcum,2));
	inds = max(0,sum(bsxfun(@gt,rnds,pcum),1))+1;
	ppnext = pgrid(inds);
	
	ppcell{jlevel+1} = reshape(ppnext,2,[]);
	level = jlevel;
	
	if doplots
	    subplot(211);
	    plot(pgrid,finitediff(pcum));
	    subplot(212);
	    mkHDPtreeplot(ppcell);

	    drawnow;
	    jj
	end
	
    end
    
    % now sample bottom level params using data (from Beta)
    pup = repmat(ppcell{nlevels}(:)',2,1);
    ppcell{nlevels+1} = max(pmin,min(pmax,betarnd(...
	alphas(nlevels+1)*pup+N1, ...
	alphas(nlevels+1)*(1-pup)+N0)));

    if doplots
	mkHDPtreeplot(ppcell);
	drawnow;
    end

    % store sample
    piSamps(:,jj) = ppcell{nlevels+1}(:);
    topsmps(jj) = ppcell{1};
    
    
end
    
    
