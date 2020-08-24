function [err,fit] = dots3DMP_fitDDM_err(param, data, options)

global call_num

ntrials = length(data.heading);

cohs = unique(data.coherence); % visual coherence levels
hdgs = unique(data.heading);

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
% a simpler, non time-dependent relationship (conf proportional to accum
% evidence) could be used

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

choice = nan(ntrials,1);
RT = nan(ntrials,1);
finalV = nan(ntrials,1);
hitBound = zeros(1,ntrials);
logOddsCorr = nan(ntrials,1);
expectedPctCorr = nan(ntrials,1);
conf = nan(ntrials,1);

modality = data.modality;
hdg = data.heading;
coh = data.coherence;
delta = data.delta;

% fundamental change: don't simulate every trial, just every unique trial
% type; simulate it K times, do IBS, then lookup-table the log-likelihood
% for each trial and sum them

% % uTrialVect = nan(ntrials,4);
% % for n = 1:ntrials
% %     uTrialVect(n,:) = [data.modality(n) data.heading(n) data.coherence(n) data.delta(n)];
% % end
% % uTrials = unique(uTrialVect,'rows');

% NO NEED FOR THIS; if continued, this would be analogous to ibs_basic. 
% instead, can just keep our sim as normal and pass this function as an
% argument to ibslike.m


for n = 1:ntrials

    switch modality(n)
        case 1
            mu = acc .* kves * sind(hdg(n)); % mean of momentary evidence
%             dv = [0, cumsum(normrnd(mu,sigmaVes,1,dur(n)-1))];
            dv = [0, cumsum(normrnd(mu,sigmaVes))];
        case 2
            mu = vel .* kvis(cohs==coh(n)) * sind(hdg(n));
%             dv = [0, cumsum(normrnd(mu,sigmaVis(cohs==coh(n)),1,dur(n)-1))];
            dv = [0, cumsum(normrnd(mu,sigmaVis(cohs==coh(n))))];
        case 3
            % positive delta defined as ves to the left, vis to the right
%             muVes = kves                      * sind(hdg(n) - delta(n)/2);
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
            dv = [0, cumsum(normrnd(mu,sigmaComb))]; % if mu is a vector
    end

    cRT = find(abs(dv)>=B, 1);
    if isempty(cRT) % did not hit bound
        RT(n) = dur(n) + Tnd;
        finalV(n) = dv(dur(n));
        hitBound(n) = 0;
    else % hit bound
        RT(n) = cRT + Tnd;
        finalV(n) = B*sign(dv(cRT));
        hitBound(n) = 1;
    end    
    choice(n) = sign(finalV(n));
    
    % use map to look up log-odds that the motion is rightward
    diffV = abs(vAxis-finalV(n));
    diffT = abs(tAxis-RT(n));
        
    thisV = find(diffV==min(diffV));
    thisT = find(diffT==min(diffT));
    logOddsCorr(n) = logOddsCorrMap(thisV(1), thisT(1));
    expectedPctCorr(n) = logistic(logOddsCorr(n)); % convert to pct corr
    conf(n) = 2*expectedPctCorr(n) - 1; % convert to 0..1

end

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

% Do ibs_basic?



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


