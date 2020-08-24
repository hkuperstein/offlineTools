%% simple monte-carlo style simulation of a 1D drift-diffusion model in the 
% context of a motion discrimination task (see e.g. Shadlen et al. 2006)
% -CF circa 2015

clear all; close all;

ntrials = 5000;

% different levels of motion strength (percent coherence)
cohs = [-0.512 -0.256 -0.128 -0.064 -0.032 0 0.032 0.064 0.128 0.256 0.512];

% delta-T (time step, in ms)
dT = 1;
% max duration of stimulus
maxdur = 2000;
timeAxis = 0:dT:maxdur;

% randomize coherence and duration for all trials (n draws with replacement)
coh = randsample(cohs,ntrials,'true')';


%%
k = 0.3; % coefficient converting coherence to units of momentary evidence
B = 25; % height of the bound, or threshold, for decision termination
sigma = 1; % standard deviation of momentary evidence


% initialize variables
dv = nan(ntrials,length(timeAxis)); % the DV as a function of n and time
choice = nan(ntrials,1); % vector to store choices: 2=right=positive, 1=left=negative
RT = nan(ntrials,1); % reaction time
finalV = nan(ntrials,1); % endpoint of the DV
hitBound = nan(ntrials,1);


%% the diffusion process
tic
for n = 1:ntrials
    % mean of momentary evidence
    mu = k * coh(n);

    % diffusion process: slow version
    if n==1 % run this once, for illustration purposes
        momentaryEvidence = nan(maxdur,1);
        for t = 1:maxdur
            momentaryEvidence(t) = randn*sigma + mu; % random sample with mean=mu, s.d.=sigma
        end
        dv(n,1) = 0;
        dv(n,2:maxdur+1) = cumsum(momentaryEvidence);
        figure; plot(dv(n,:)); hold on; 
        plot(1:length(dv),ones(1,length(dv))*B,'g-');
        plot(1:length(dv),ones(1,length(dv))*-B,'r-');
        xlabel('Time (ms)'); ylabel('Accum. evidence (DV)');
        % (evidence is shown continuing to accumulate
        % past the bound, although it realy stops there)
    end    

    % faster version: does not require a FOR loop over the variable t
    dv(n,:) = [0, cumsum(normrnd(mu,sigma,1,maxdur))];
        
    cRT = find(abs(dv(n,:))>=B, 1);
    if isempty(cRT) % did not hit bound
        RT(n) = maxdur;
        finalV(n) = dv(RT(n));
        hitBound(n) = 0;
    else % hit bound
        RT(n) = cRT;
        finalV(n) = B*sign(dv(n,RT(n)));
        hitBound(n) = 1;
    end    
    choice(n) = sign(finalV(n));  
end
toc



pCorrect_total = (sum(choice==1 & coh>0) + sum(choice==-1 & coh<0)) / ntrials;




%% plot proportion "rightward" (choice=1) and reaction time as a function of coherence

for c = 1:length(cohs)
    I = coh==cohs(c);
    pRight(c,1) = sum(I & choice==1) / sum(I);
    meanRT(c,1) = mean(RT(I));
end

figure; plot(cohs,pRight(:,1),'bo-');
xlabel('Motion strength (%coh)'); ylabel('Proportion rightward choices');

figure; plot(cohs,meanRT(:,1),'go-');
xlabel('Motion strength (%coh)'); ylabel('Reaction time (ms)');




