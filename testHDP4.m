% testHDP4.m
%
% Examine performance of Gibbs sampler.

addpath HDPtools/  % add path 

nlevels = 3;
alphas = 100;%20*1.5.^(0:nlevels); % concentration parameters for each level
pmu = .5;  % mean spike probability
nstates = 2^nlevels;
nsamps = 5;  % number of samples to draw

%% Create Markov model
%nstates = 8*nstates;
% Set true transition probabilities by hand
ps = (0.5:(-0.4/(nstates-1)):0.1)';
%ps = [.6, .61, .5,.5,.2,.2, .25,.23]';

% Insert into transition matrix
[T,Tpairs,ppairs] = mkBinaryMCtransitionMatrix(ps);

subplot(222);
imagesc(1:nstates,1:nstates,T); colormap gray; axis image;
title('Transition matrix');
statelabels = fliplr(dec2bin(0:nstates-1));
xtcks = unique(round(linspace(0,nstates,5)));
set(gca,'ytick',1:nstates, 'xtick',xtcks(2:end));
set(gca,'yticklabel',statelabels,'xaxislocation','top',...
    'tickdir','out','xticklabel',statelabels(xtcks(2:end),:));
box off;
xlabel('state t+1');
ylabel('state t');

% Compute stationary distribution
u = compStationaryDistFromT(T);

%% Simulate from Markov chain

nsamps = 1000;
X = zeros(nsamps,1);
X(1) = 1;  % initialize in first state (for now).
for j = 2:nsamps
    X(j) = randsample(Tpairs(X(j-1),:),1,true,ppairs(X(j-1),:));
end
spktrain = zeros(nsamps,1);
spktrain(X(1:end)>nstates/2) = 1;

subplot(426);
plot(1:nsamps,X,'.-');
set(gca,'ylim',[1 nstates],'ytick',1:nstates,'yticklabel',statelabels);
box off; set(gca,'tickdir','out');
ylabel('state');

subplot(428);
plot(1:nsamps,spktrain,'.-');
box off; set(gca,'tickdir','out','ytick',[0 1]);
set(gcf,'color','w');
ylabel('spike');


%% Count empirical state transition and marginal state distributions
[Tcounts,Tempir,pstate_empir] = countMCstateTransitions(X,nlevels);
pasymp_empir = compStationaryDistFromT(Tempir);  % eigenvector estimate of stationary distribution
pshat = unmkBinaryMCtransitionMatrix(Tempir); % estimated transition probabilities

subplot(224);
imagesc(1:nstates,1:nstates,Tempir); axis image;
set(gca,'ytick',1:nstates, 'xtick',xtcks(2:end));
set(gca,'yticklabel',statelabels,'xaxislocation','top',...
    'tickdir','out','xticklabel',statelabels(xtcks(2:end),:));
box off;
xlabel('state t+1');
ylabel('state t');
title('Empirical transition matrix');


%% Run Gibbs sampler
npgrid = 100;
ngibbssamps =1000;
doplots = 0;  % boolean for whether or not to visualize Gibbs samples
tic;
[psmps,pusmps] = gibbsSampleBinaryHDP(pmu,alphas,Tcounts,ngibbssamps,npgrid,doplots);
toc;
ps_bls = mean(psmps')'; % BLS estimate for p
ps_eb = std(psmps')';  % 1SD error bars 

% Compute asymptotic distribution associated with BLS estimate
[Tbls] = mkBinaryMCtransitionMatrix(ps_bls);
ubls = compStationaryDistFromT(Tbls);

%% Examine Gibbs samps

subplot(425);
plot(1:nstates,ps,'.-k', 1:nstates,pshat,'o-b', 1:nstates, ps_bls,'o-r');
hold on;
errorbar(1:nstates,ps_bls,2*ps_eb,'r');
hold off;
set(gca,'ylim', [0 max([ps;pshat])*1.2]);
ylabel('P(1|contexts)');
legend('true','plugin','BLS');
set(gca,'xlim',[1 nstates]);
box off;

Errs_TransitionPs = [norm(ps-pshat), norm(ps-ps_bls)]

subplot(427);
plot(1:nstates,u,'.-k',1:nstates,pasymp_empir,'o-',1:nstates,ubls, 'ro-');
ylabel('marginal: P(states)');
set(gca,'xlim',[1 nstates]);
xlabel('state');
box off;

Errs_StationaryPs = [norm(u-pasymp_empir), norm(u-ubls)]


%% Examine Gibbs samples
subplot(221);
plot(1:ngibbssamps,psmps',1:ngibbssamps,ones(ngibbssamps,1)*ps', 'k--');
title('Gibbs samples');
ylabel('p(1|state)');
xlabel('sample #');


%% Look at entropy rate estimates

mu = mean(spktrain);
Hub = -mu*log2(mu)-(1-mu)*log2(1-mu);
Hplug = compBinaryMCentropyrate(pshat,pstate_empir);
Hblsplug = compBinaryMCentropyrate(ps_bls,ubls);
Htrue = compBinaryMCentropyrate(ps,u);

Hestims = [Htrue,Hub,Hplug, Hblsplug]
