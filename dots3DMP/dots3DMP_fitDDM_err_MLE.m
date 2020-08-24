function [err,fit] = dots3DMP_fitDDM_err_MLE(param, data, options)

global call_num

ntrials = length(data.heading);

nsimtrl = 100; %sample size of Monte Carlo simulation per trial

cohs = unique(data.coherence); % visual coherence levels
hdgs = [-10 -5 -2.5 -1.25 -eps eps 1.25 2.5 5 10]; % heading angles
                            % (map fn seems to require even number of diff 
                            % levels (or just no zero), so we use +/- eps)
duration = 2000; % stimulus duration (ms)
dur = ones(ntrials,1) * duration;

kves = param(1);
kvis = param(2)*cohs;
B = abs(param(3)); % don't accept negative bound heights

% what about sigma? assume 1 for now (later may need not just separate 
% sigmas but full Drugowitsch paramaterization):
sigmaVes = 1; % std of momentary evidence
sigmaVis = [1 1]; % allow for separate sigmas for condition, coherence

Tnd = 300; % non-decision time (or make this the mean of a dist?)

maxdur = duration;
% assume the mapping is based on an equal amount of experience with the 
% *three* levels of reliability (ves, vis-low, vis-high) hence k and sigma
% are their averages
k = mean([kves kvis']);
sigma = mean([sigmaVes sigmaVis]);
[~, ~, logOddsCorrMap, tAxis, vAxis] = makeLogOddsCorrMap_3DMP(hdgs,k,B,sigma,maxdur,0);
% uses Fokker-Planck equation to propagate the probability density of the DV,
% as in Kiani & Shadlen 2009. Required for readout of confidence, although
% a simpler heuristic could be used (conf proportional to accum evidence)

% create acceleration and velocity profiles (arbitrary for now)
% SJ 04/2020
% Hou et al. 2019, peak vel = 0.37m/s, SD = 210ms
vel = normpdf(1:duration,duration/2,210);
vel = 0.37*vel./max(vel);
acc = gradient(vel)*1000; % multiply by 1000 to get from m/s/ms to m/s/s

% normalize
vel = vel./max(vel);
acc = abs(acc./max(acc)); % (and abs)


%% bounded evidence accumulation

choice           = nan(ntrials,nsimtrl);
RT               = nan(ntrials,nsimtrl);
finalV           = nan(ntrials,nsimtrl);
hitBound         = false(ntrials,nsimtrl);
logOddsCorr      = nan(ntrials,nsimtrl);
expectedPctCorr  = nan(ntrials,nsimtrl);
conf             = nan(ntrials,nsimtrl);

hdg = data.heading;
coh = data.coherence;
delta = data.delta;
modality = data.modality;

tic
for n = 1:ntrials

    switch modality(n)
        case 1
            mu = acc .* kves * sind(hdg(n)); % mean of momentary evidence
%             dv = [0, cumsum(normrnd(mu,sigmaVes,1,dur(n)-1))];
%             dv = [0, cumsum(normrnd(mu,sigmaVes))];

            dv = [zeros(nsimtrl,1) cumsum(repmat(mu,nsimtrl,1) + randn(nsimtrl,maxdur)*sigmaVes, 2)];

        case 2
            mu = vel .* kvis(cohs==coh(n)) * sind(hdg(n));
%             dv = [0, cumsum(normrnd(mu,sigmaVis(cohs==coh(n)),1,dur(n)-1))];
%             dv = [0, cumsum(normrnd(mu,sigmaVis(cohs==coh(n))))];
            dv = [zeros(nsimtrl,1) cumsum(repmat(mu,nsimtrl,1) + randn(nsimtrl,maxdur)*sigmaVis(cohs==coh(n)), 2)];

        case 3
            % positive delta defined as ves to the left, vis to the right
            muVes = acc .* kves               * sind(hdg(n) - delta(n)/2);
            muVis = vel .* kvis(cohs==coh(n)) * sind(hdg(n) + delta(n)/2);
            
            wVes = sqrt( kves^2 / (kvis(cohs==coh(n))^2 + kves^2) );
            wVis = sqrt( kvis(cohs==coh(n))^2 / (kvis(cohs==coh(n))^2 + kves^2) );

            mu = wVes.*muVes + wVis.*muVis;
            % here the DV is a sample from a dist with mean = weighted sum
            % of means. thus the variance is the weighted sum of variances
            % (error propagation formula):
            sigmaComb = sqrt(wVes.^2 .* sigmaVes^2 + wVis.^2 .* sigmaVis(cohs==coh(n))^2); % assume zero covariance
%             dv = [0, cumsum(normrnd(mu,sigmaComb,1,dur(n)-1))];
%             dv = [0, cumsum(normrnd(mu,sigmaComb))]; % if mu is a vector
            dv = [zeros(nsimtrl,1) cumsum(repmat(mu,nsimtrl,1) + randn(nsimtrl,maxdur)*sigmaComb, 2)];

    end

    % get RT, choice, and confidence from DV
    [~,cRT] = max(abs(dv)>=B,[],2);
    hitBound(n,:) = cRT>1;
    cRT(cRT==1) = dur(n); % did not hit bound, cRT==dur
    
    RT(n,:) = cRT + Tnd;
    
    RTind  = sub2ind(size(dv),(1:nsimtrl)',cRT);

    finalV(n,~hitBound(n,:)) = dv(RTind(~hitBound(n,:)));
    finalV(n,hitBound(n,:))  = B.*sign(dv(RTind(hitBound(n,:))));
%     choice(n,:) = sign(finalV(n,:));
    
    % use map to look up log-odds that the motion is rightward
    [~,thisV] = min(abs(bsxfun(@minus,vAxis,finalV(n,:)')),[],2);
    [~,thisT] = min(abs(bsxfun(@minus,tAxis',cRT)),[],2);
    
    logmapinds = sub2ind(size(logOddsCorrMap),thisV,thisT);
    
    logOddsCorr(n,:) = logOddsCorrMap(logmapinds);
%     expectedPctCorr(n,:) = logistic(logOddsCorr(n,:));
%     conf(n,:) = 2*expectedPctCorr - 1;

end
toc
keyboard

expectedPctCorr = logistic(logOddsCorr);
conf = 2*expectedPctCorr - 1;
    

choice = sign(finalV);
choice(choice==0) = sign(randn); % not needed under usual circumstances
choice(choice==1) = 2; choice(choice==-1) = 1; % 1=left, 2=right

% output var 'fit' gets the same values as data for the conditions, but the
% simulated trial outcomes for the observables!
fit.heading = data.heading;
fit.coherence = data.coherence;
fit.modality = data.modality;
fit.delta = data.delta;
fit.choice = choice;
fit.RT = RT;
fit.conf = conf;


%% calculate error

% What we want is to use maximum-likelihood estimation (MLE), so we need a
% way to calculate the likelihood.
% Here's a refresher (see also likelihood_tutorial.m):

% The likelihood is often written P(data|params), which looks like a 
% conditional density where data varies and params is fixed, but it's
% really kind of the opposite: it's a function of the params for fixed data.

% Confusingly, we don't actually need an equation for the full likelihood
% *function* -- that's the Nparams-dimensional surface which is mapped out
% over many iterations of our fitting routine. The optimization just needs
% to sample the surface sufficiently well to find the maximum.

% What each iteration needs is a functional form for the conditional
% density p(data|params), where data is a set of N trials x 3 observations,
% choice, RT, and confidence. In other words we need a concise way to
% describe the stochastic process that gives rise to the observations,
% given a set of params.

% In likelihood_tutorial.m, the stochastic process is Poisson, the
% observation is a spike count, and params is one-dimensional: some
% stimulus value theta. The probability p(obs|theta) for a given theta
% is either calculated from the equation for the Poisson probability
% mass function, or estimated using a 'histogram method' by generating many
% samples from this distribution.

% This approach seems like it might be amenable to what we're doing here,
% since we're generating the model fits/predictions through monte carlo
% simulation anyway. But there's a key difference: In likelihood_tutorial,
% there was only a scalar observation (the spike count) and an IID set of
% trials (poissrnd(theta)). Here, we have a vector of observations (choice,
% RT, conf) and a non-iid set of trials, each defined by the 4-vector:
% [heading, coh, delta, modality].

% In the Kiani DDM fitting code, it was possible to accumulate (log)
% likelihood trial by trial, because for any combination of coh and
% duration, the model specified a probability of rightward, leftward, and
% sure bet *for that trial*. To replicate this in our monte carlo sim-style
% fitting, we'd have to simulate every combination of hdg, mod, coh, delta
% many times to estimate these probabilities. This is probably too slow to
% be viable, but I may try it anyway.

% normalize before calculating error:
% shift choice back to 0..1, normalize RT, and conf already is 0..1
choiceM = choice - 1;
choiceD = data.choice - 1;

RTm = RT/max(RT);
RTd = data.RT/max(data.RT);

% err_choice = sum((choiceM-choiceD).^2);
% err_RT = sum((RTm-RTd).^2);
% err_conf = sum((conf-data.conf).^2);


% calculate expected probability of rightward choice given model params
% Eq. from Shadlen et al. 2006 book chapter:
% Pr_model = 1 ./ (1 + exp(-2*k*B*stimval));
% Pr_model = 1 ./ (1 + exp(-2*B*meanMu)); 

% 'stimval' needs to be some function of k and mu (which varies over time,
% unlike RDM task version)

% estimate Pright from Monte Carlo simulation
Pr_model = sum(choiceM,2)/nsimtrl;

% kluge to avoid taking log(0) below
Pr_model(Pr_model==0) = min(Pr_model(Pr_model~=0));
Pr_model(Pr_model==1) = max(Pr_model(Pr_model~=1));

choiceD = logical(choiceD);
LL_choice = sum(log(Pr_model(choiceD))) + sum(log(1-Pr_model(~choiceD)));

err_choice = -LL_choice; % we want to minimize NLL
err = err_choice; % placeholder

% TODO confidence and RT error calculations
% a few options here... 
% 1. fit distributions
% 2. fit quantiles
% 3. fit mean (for RT only?)

% total err should then be
% err = err_choice + err_RT + err_conf;


%% print progress report!
fprintf('\n\n\n****************************************\n');
fprintf('run %d\n', call_num);
fprintf('\tkves= %g\n\tkvisMult= %g\n\tB= %g\n', param(1), param(2), param(3));
fprintf('err: %f\n', err);
if options.ploterr && strcmp(options.fitMethod,'fms')
    figure(options.fh);
    plot(call_num, err, '.','MarkerSize',12);
    drawnow;
end
call_num = call_num + 1;



end


