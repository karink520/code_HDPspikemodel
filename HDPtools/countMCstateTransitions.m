function [Tcounts,T_empir,pstate_empir] = countMCstateTransitions(X,depth)
% [Tcounts,T_empir,pstate_empir] = countMCstateTransitions(X,depth)
%
% Count Markov Chain state transitions, and estimate state transition
% matrix T and marginal state probabilities

Nstates = 2^depth; % number of states
Nh = Nstates/2;    % half the number of states


% Compute unique symbols for each pair of states 
jointlabels = X(1:end-1)-1 + (2^depth)*(X(2:end)-1);  % map state pairs to unique symbols

% Now do the same thing for matrix sparsity pattern
I = [(1:Nstates)'; (1:Nstates)'];
J = reshape([1:Nh;1:Nh],[],1); J = [J;J+Nh];
Tlabels = (I-1) + (2^depth)*(J-1);

% Now count state transitions
Tcounts = hist(jointlabels,Tlabels);


if nargout > 1
    % Estimate marginal state probabilities (excluding last state)
    pstate_empir = hist(X(1:end-1),1:Nstates)';
    
    % Insert counts into empirical T matrix estimate
    T_empir = sparse(I,J,Tcounts);
    
    % Normalize by # of observations to obtain plugin T esimate
    ii = find(pstate_empir);
    T_empir(ii,:) = bsxfun(@rdivide,T_empir(ii,:),pstate_empir(ii));
    
    pstate_empir = pstate_empir/sum(pstate_empir); % normalize to 1
end


Tcounts = reshape(Tcounts,[],2);
