%% monte-carlo simulation of Harry's bandit-dots task:
% a two armed bandit where an bandit choice is followed by an RDM decision
% with confidence and RT. The reward is drawn from one of three overlapping
% distributions: low, medium, and high. Low when both the bandit and dots
% choices are wrong, high when both are correct, and medium when either is
% correct and the other incorrect. The correct bandit choice is random at
% the beginning of the block of ~40-60 trials, and changes at some uncued
% time near the middle of the block. In addition, expected values of the
% options can be modulated by changing, say, the prior over (signed)
% coherence, in addition to mean and variance of reward distributions.


% Questions we want to address with this:
% how does confidence influence explore/exploit decisions? (irrationally?)
% how does perceptual uncertainty affect the learning rate?

% Harry: "This proposal aims to investigate interactions between value and
% perceptual uncertainties, how decision confidence is affected by both,
% and how confidence may mediate exploration and exploitation in multi-
% option choices.

% Subjects begin by making a choice between two options: A and B. The
% options have different underlying reward distributions with unequal means
% and equal variances. To receive a reward drawn from the chosen option?s
% distribution, subjects are required to correctly discriminate the
% direction of motion in a noisy visual stimulus, displayed immediately
% after the A vs. B choice. They will also report their confidence in
% their perceptual choice at this stage. Subjects will complete 30 trials
% per block, for a total of 10 experimental blocks. Their goal is to
% maximize their reward, which will correlate with a monetary reward,
% intended to motivate subjects to optimize. The PURA funding will
% assist in this regard.

% Data Analysis: On each trial we will record the subject's choices,
% reaction times, and confidence in the perceptual decision. These data
% will be used to develop and validate a combined drift-diffusion and
% reinforcement learning model, including the extraction of value learning
% rates using reaction times (Ballard and McClure, 2019), a novel method
% which this project would be the first in the field to demonstrate.
% Comparing choice progression and confidence reports will then elucidate
% the role of confidence in both exploratory behavior and changes in
% learning rate."

% So the idea is RT on the A/B choice is related to value. It should also
% be affected by confidence on the previous trial(s). It may do so in an
% irrational way, like the more confident you are the faster you go and 
% thus the greater your value estimate is for the chosen option -- even 
% independently of the reward you received on that high-confidence trial

% What I don't have is a hypothesis. All these things may be related,
% conf, RT, learning rate -- but it's bound to be complex and hard to
% constrain with data. We need a specific hypothesis we can test.

% A preliminary question I have is, if not softmax, what is the decision
% rule for the initial (bandit) choice? We need a DV that implements a
% dividend for exploration; has that been done?
% Closest thing may be Purcell, where an across-trial DV accumulates up to
% a 'switch' decision. I'd rather have it be a graded or at least noisy
% thing, maybe sampled from a dist with mean that increments or decrements
% every trial but is still noisy and accumulates...

% Let's think about whether it is possible to ask the deep question
% underlying the above thought: do brains generate random numbers to
% facilitate exploration, or not?

% The two accounts would predict different things given a manipulation of
% noise. if noise accounts for exploration, then reducing noise should
% reduce it. If intentional randomization accounts for it, then it won't.
% How do you manipulate noise? Dstoch won't work, because we're talking
% about noise in the DV that governs the A/B choice. And whatever we do
% must permit a distinction between noise and the value itself. Difficult
% to see how confidence can serve this role.

% What about casting this as information seeking? you always choose the
% item with greater value, you don't flip a weighted coin, but 'value'
% includes a term for the value of information, which will naturally
% promote explore choices.



%% simulate bandit choice


%% option 1: Q learning / Rescorla-Wagner style (action values + softmax)

% Findling et al (Wyart)
% Q(t) = Q(t-1) + alpha*(r(t-1) - Q(t-1)) + epsilon(t) % expected value of chosen action a(t-1)
% a(t) = Bernoulli(1 / 1 + exp(-beta(Q(t,A)-Q(t,B)) - xi*sign(a(t)-1))) % softmax with choice histeresis

% OR

% Bari et al (Cohen)
% if c(t) = l,
% Ql(t+1) = zeta*Ql(t) + alpha*(R(t)-Ql(t))
% Qr(t+1) = zeta*Qr(t)
% if c(t) = r,
% Ql(t+1) = zeta*Ql(t)
% Qr(t+1) = zeta*Qr(t) + alpha*(R(t)-Qr(t))    
% (i.e. only the chosen option is updated by alpha*delta)

% p(c(t)=r) = 1 / 1+exp(-beta(Qr(t)-Ql(t)+b+kappa*a(t-1))) % basically the same as above
% p(c(t)=l) = 1 - p(c(t)=r)

% initialy it seemed like we won't need the Bari formulation because
% expected value of the two options should not vary independently: when one
% is correct, the other is incorrect. However, this may not be true for
% actual subjects, not to mention that we may want to allow for things like
% perceived accuracy/confidence, distribution of cohs, etc to differ
% between the options.
% the other way to think about it is that the value representations are not
% memoryless. any evidence you have that A is correct should persist even
%if you switch to B for a while. this is why a simple switch evidence model
%(with hard resetting) like Purcell doesn't seem quite right either...


% AND after Bartolo @DMJC, don't we need to consider model-based (Bayesian)
% and compare to model-free?


clear all; close all;
plotflag = 1;

ntrials = 200; % check w harry

% initialize bandit params
zeta = 1; % 1 minus forgetting rate
alpha1 = 1; % learning rate (starting value * 1/average conf; see below)
beta = 0.8; % inverse temperature
b = 0; % bias
kappa = 0; % histeresis

Ql = nan(ntrials,1); % expected value of left
Qr = nan(ntrials,1); % expected value of right
Pr = nan(ntrials,1); % probability of rightward choice (Pl = 1-Pr)
c = nan(ntrials,1); % bandit choice: -1 for left, +1 for right
s = nan(ntrials,1); % 'state': which option is correct?
R = nan(ntrials,1); % actual reward delivered
alpha = nan(ntrials,1); % learning rate (will vary trial to trial e.g. based on confidence) 

% reward distributions
rewMu = [5 10 15];
rewSigma = [3 3 3];

if plotflag==2
    x = 0:0.1:20;
    figure; plot(x,normpdf(x,rewMu(1),rewSigma(1)),'r',x,normpdf(x,rewMu(2),rewSigma(2)),'b',x,normpdf(x,rewMu(3),rewSigma(3)),'g');
%         Q for harry: was dist 3 truncated at the high end??
%         what was the SD?
end


% dots task params
cohs = [-0.512 -0.256 -0.128 -0.064 -0.032 0 0.032 0.064 0.128 0.256 0.512];
coh = randsample(cohs,ntrials,'true')';
dT = 1; % ms
maxdur = 2000; % max viewing duration (RT task)
timeAxis = 0:dT:maxdur;
k = 0.25; % drift rate coeff (conversion from %coh to units of DV)
B = 25; % bound height
mu = k*coh; % mean of momentary evidence (drift rate)
sigma = 1; % SD of momentary evidence
theta = 1; % threshold for Psure choice in units of log odds correct (Kiani & Shadlen 2009)

% Kiani map [need to get 2014 version!]
[logOddsMapR, logOddsMapL, logOddsCorrMap, tAxis, vAxis] = makeLogOddsCorrMap_smooth(k,B,sigma,theta,timeAxis,0);

% non-decision time (truncated normal dist)
TndMean = 300;
TndSD = 50;
TndMin = TndMean/2;
TndMax = TndMean+TndMin;
Tnd = zeros(ntrials,1);
for t = 1:ntrials
    while Tnd(t)<=TndMin || Tnd(t)>=TndMax % simple trick for truncating
        Tnd(t) = round(normrnd(TndMean,TndSD));
    end
end


%%

% initialize bandit
Ql(1) = 10; Qr(1) = 10; R(1) = 10; % mean of middle dist
c(1) = sign(rand); % first real trial is actually t=2, but needs a starting point 
alpha(2) = alpha1*0.7; % alpha(1) stays undefined
s(:) = 1; % temp; will be sign(rand)
s(round(ntrials/2):end) = -s(1); % temp; will be a finite switching prob in middle trials with some refractory period

% initialize dots
dotsChoice = nan(ntrials,1); % choices (left = -1, right = 1);
DT = nan(ntrials,1); % decision time (or time-to-bound for fixed/variable duration task)
RT = nan(ntrials,1); % measured reaction time, which includes non-decision time
expectedPctCorr = nan(ntrials,1);
conf = nan(ntrials,1);
hitBound = nan(ntrials,1);
logOddsCorr = nan(ntrials,1);
corBandit = nan(ntrials,1);
corDots = nan(ntrials,1);

for t = 2:ntrials
    
    % action values from Q learning (Bari et al)
    Ql(t) = zeta * Ql(t-1) + (c(t-1)==-1) * alpha(t) * (R(t-1) - Ql(t-1));
    Qr(t) = zeta * Qr(t-1) + (c(t-1)==1)  * alpha(t) * (R(t-1) - Qr(t-1));
        % by convention, let alpha(t) be the learning rate for the given
        % trial, even though it is based on confidence on the previous
        % trial and modulates effect of RPE on previous trial
    Pr(t) = logistic(beta * (Qr(t) - Ql(t)) + b + kappa*c(t-1)); % softmax
    c(t) = (rand>(1-Pr(t))) * 2 - 1; % weighted coin flip
    corBandit(t) = c(t)==s(t); % correct if choice matches state
           
    % now the dots task
    
    % the diffusion process
    dv = [0, cumsum(normrnd(mu(t),sigma,1,maxdur))];

    cRT = find(abs(dv)>=B, 1, 'first');
    if isempty(cRT)
        RT(t) = maxdur;
        dvEnd = dv(RT(t));
        hitBound(t) = 0;
    else
        RT(t) = cRT;
        dvEnd = B*sign(dv(RT(t))); % cap bound crossings at bound value
        hitBound(t) = 1;
    end
   
    % choice
    if dvEnd==0
        dotsChoice(t) = sign(rand-0.5);
    else
        dotsChoice(t) = sign(dvEnd);
    end
    corDots(t) = sign(dotsChoice(t))==sign(coh(t));
    dotsChoice(t) = (dotsChoice(t)+1)/2; % convert to 0/1
    
    % look up the expected log odds corr for this trial's V and T
    whichV = find(abs(vAxis-dvEnd)==min(abs(vAxis-dvEnd)));
    whichT = find(abs(tAxis-RT(t))==min(abs(tAxis-RT(t))));   
    logOddsCorr(t) = logOddsCorrMap(whichV(1), whichT(1));

    % confidence rating
    expectedPctCorr(t) = logistic(logOddsCorr(t)); % convert to pct corr
    conf(t) = 2*expectedPctCorr(t) - 1; % convert to 0..1    
    
    % reward
    rewdist = 1+corDots(t)+corBandit(t);
    R(t) = round(normrnd(rewMu(rewdist),rewSigma(rewdist))); % truncate?
    
    % adjust alpha for next trial
    if t<ntrials
        alpha(t+1) = conf(t)*alpha1;
    end
    
end

figure; plot(c); ylim([-1.25 1.25]);

%% plot proportion "rightward" (choice=1) and reaction time as a function of coherence

for C = 1:length(cohs)
    I = coh==cohs(C);
    pRight(C,1) = sum(I & dotsChoice==1) / sum(I);
    meanRT(C,1) = mean(RT(I));
    meanConf(C,1) = mean(conf(I));
end

figure(122); clf; set(gcf,'Position',[900   465   480   750]);
subplot(3,1,1); plot(cohs,pRight(:,1),'co-'); ylim([0 1]);
xlabel('Motion strength (%coh)'); ylabel('Proportion rightward choices');

subplot(3,1,2); plot(cohs,meanRT(:,1),'go-');
xlabel('Motion strength (%coh)'); ylabel('Reaction time (ms)');

subplot(3,1,3); plot(cohs,meanConf(:,1),'ro-'); ylim([0 1]);
xlabel('Motion strength (%coh)'); ylabel('Confidence rating');


%% option 2: accumulation of switch evidence (Sarafyazd & Jazayeri, Purcell+Kiani)

clear

ntrials = 200; % check w harry
% reward distributions
rewMu = [5 10 15];
rewSigma = [3 3 3];

% initialize switching decision vairable
ve = nan(ntrials,1);
ve(1)=0;

%initialize boundary parameters
sigma_e = 1;  % switching noise
lambda = 0; % leakage for past evidence
q=-inf; % negative switch evidence
B0 = 0.1; % switching parameters
gamma = 2; % switching parameters
omega = 1; % switching parameters
H = nan(ntrials,1);
Be=3; % switching bound

%initialize switching evidence
mue = nan(ntrials,1); %switch evidence
mue(1)=0;

%initialize other parameters
fb_p= zeros(ntrials,1); %positive feedback
choice_en =  zeros(ntrials,1); %choice of the environment
choice_en(1) = sign(rand);
R = nan(ntrials,1); % actual reward delivered
RT = nan(ntrials,1); % measured reaction time, which includes non-decision time
hitBound = nan(ntrials,1);
corDots = nan(ntrials,1);
corBandit = nan(ntrials,1);
dotsChoice = nan(ntrials,1);
s = nan(ntrials,1); % 'state': which option is correct?
s(:) = 1; % temp; will be sign(rand)
s(round(ntrials/2):end) = -s(1); % temp; will be a finite switching prob in middle trials with some refractory period

% dots task params
cohs = [-0.512 -0.256 -0.128 -0.064 -0.032 0 0.032 0.064 0.128 0.256 0.512];
coh = randsample(cohs,ntrials,'true')';
dT = 1; % ms
maxdur = 2000; % max viewing duration (RT task)
timeAxis = 0:dT:maxdur;
k = 0.25; % drift rate coeff (conversion from %coh to units of DV)
B = 25; % bound height
mu = k*coh; % mean of momentary evidence (drift rate)
sigma = 1; % SD of momentary evidence
theta = 1; % threshold for Psure choice in units of log odds correct (Kiani & Shadlen 2009)

% Kiani map [need to get 2014 version!]
[logOddsMapR, logOddsMapL, logOddsCorrMap, tAxis, vAxis] = makeLogOddsCorrMap_smooth(k,B,sigma,theta,timeAxis,0);


for i = 2: ntrials
    % now the dots task
    % the diffusion process
    dv = [0, cumsum(normrnd(mu(i),sigma,1,maxdur))];

    cRT = find(abs(dv)>=B, 1, 'first');
    if isempty(cRT)
        RT(i) = maxdur;
        dvEnd = dv(RT(i));
        hitBound(i) = 0;
    else
        RT(i) = cRT;
        dvEnd = B*sign(dv(RT(i))); % cap bound crossings at bound value
        hitBound(i) = 1;
    end
   
    % choice
    if dvEnd==0
        dotsChoice(i) = sign(rand-0.5);
    else
        dotsChoice(i) = sign(dvEnd);
    end
    corDots(i) = sign(dotsChoice(i))==sign(coh(i));
    dotsChoice(i) = (dotsChoice(i)+1)/2; % convert to 0/1
    
    whichV = find(abs(vAxis-dvEnd)==min(abs(vAxis-dvEnd)));
    whichT = find(abs(tAxis-RT(i))==min(abs(tAxis-RT(i))));   
    logOddsCorr(i) = logOddsCorrMap(whichV(1), whichT(1));

    % confidence rating
    expectedPctCorr(i) = logistic(logOddsCorr(i)); % convert to pct corr
    mue(i) = log(1./(1-expectedPctCorr(i))); %switch evidence

    if ve(i-1) > Be
        ve(i)=0;
        choice_en(i) = -choice_en(i-1);
        corBandit(i) = choice_en(i)==s(i);
        rewdist = 1+corDots(i)+corBandit(i);
        R(i) = round(normrnd(rewMu(rewdist),rewSigma(rewdist))); % truncate? 
    
    else
        choice_en(i) = choice_en(i-1);
        corBandit(i) = choice_en(i)==s(i); % correct if choice matches state
        rewdist = 1+corDots(i)+corBandit(i);
        R(i) = round(normrnd(rewMu(rewdist),rewSigma(rewdist))); % truncate? 
        
        p = [normpdf(R(i),rewMu(1),rewSigma(1));normpdf(R(i),rewMu(2),rewSigma(2));normpdf(R(i),rewMu(3),rewSigma(3))];
        max_p = find(p==max(p));
        if max_p==3
            mue(i) = q;
            fb_p(i) = 1;
        end
        ve(i) = ve(i-1)+ normrnd((mue(i)-lambda*ve(i-1)),sigma_e);
        %h = p- (1-expectedPctCorr(i));

    end
    if ve(i)<0
        ve(i)=0;
    end
end


figure
plot(ve, '*-')
%hold on
%plot(fb_p, 'o')
hold on
plot(choice_en, 'o')
hold on 
plot(s, 'k-')
hold on
plot(1:1:ntrials, repmat(Be, ntrials, 1), '-')
legend('switching variable', 'environment choice', 'environment', 'switching bound')



%% option 3: accumulation of switch evidence with hazard rate (Sarafyazd & Jazayeri, Purcell+Kiani)
clear

ntrials = 200; % check w harry
% reward distributions
rewMu = [5 10 15];
rewSigma = [3 3 3];

% initialize switching decision vairable
ve = nan(ntrials,1);
ve(1)=0;

%initialize boundary parameters
sigma_e = 1;  % switching noise
lambda = 0; % leakage for past evidence
q=-inf; % negative switch evidence
B0 = 0.1; % switching parameters
gamma = 2; % switching parameters
omega = 1; % switching parameters
H = nan(ntrials,1);
Be=3; % switching bound
Be0 = Be;

%initialize switching evidence
mue = nan(ntrials,1); %switch evidence
mue(1)=0;

%initialize other parameters
fb_p= zeros(ntrials,1); %positive feedback
choice_en =  zeros(ntrials,1); %choice of the environment
choice_en(1) = sign(rand);
R = nan(ntrials,1); % actual reward delivered
RT = nan(ntrials,1); % measured reaction time, which includes non-decision time
hitBound = nan(ntrials,1);
corDots = nan(ntrials,1);
corBandit = nan(ntrials,1);
dotsChoice = nan(ntrials,1);
s = nan(ntrials,1); % 'state': which option is correct?
s(:) = 1; % temp; will be sign(rand)
s(round(ntrials/2):end) = -s(1); % temp; will be a finite switching prob in middle trials with some refractory period

% dots task params
cohs = [-0.512 -0.256 -0.128 -0.064 -0.032 0 0.032 0.064 0.128 0.256 0.512];
coh = randsample(cohs,ntrials,'true')';
dT = 1; % ms
maxdur = 2000; % max viewing duration (RT task)
timeAxis = 0:dT:maxdur;
k = 0.25; % drift rate coeff (conversion from %coh to units of DV)
B = 25; % bound height
mu = k*coh; % mean of momentary evidence (drift rate)
sigma = 1; % SD of momentary evidence
theta = 1; % threshold for Psure choice in units of log odds correct (Kiani & Shadlen 2009)

% Kiani map [need to get 2014 version!]
[logOddsMapR, logOddsMapL, logOddsCorrMap, tAxis, vAxis] = makeLogOddsCorrMap_smooth(k,B,sigma,theta,timeAxis,0);

count_error = 0;
h_subjective = [];
hazard_temp = 0;

for i = 2: ntrials
    % now the dots task
    % the diffusion process
    dv = [0, cumsum(normrnd(mu(i),sigma,1,maxdur))];

    cRT = find(abs(dv)>=B, 1, 'first');
    if isempty(cRT)
        RT(i) = maxdur;
        dvEnd = dv(RT(i));
        hitBound(i) = 0;
    else
        RT(i) = cRT;
        dvEnd = B*sign(dv(RT(i))); % cap bound crossings at bound value
        hitBound(i) = 1;
    end
   
    % choice
    if dvEnd==0
        dotsChoice(i) = sign(rand-0.5);
    else
        dotsChoice(i) = sign(dvEnd);
    end
    corDots(i) = sign(dotsChoice(i))==sign(coh(i));
    dotsChoice(i) = (dotsChoice(i)+1)/2; % convert to 0/1
    
    whichV = find(abs(vAxis-dvEnd)==min(abs(vAxis-dvEnd)));
    whichT = find(abs(tAxis-RT(i))==min(abs(tAxis-RT(i))));   
    logOddsCorr(i) = logOddsCorrMap(whichV(1), whichT(1));

    % confidence rating
    expectedPctCorr(i) = logistic(logOddsCorr(i)); % convert to pct corr
    mue(i) = log(1./(1-expectedPctCorr(i))); %switch evidence

    if ve(i-1) > Be(end)
        ve(i)=0;
        choice_en(i) = -choice_en(i-1);
    else
        choice_en(i) = choice_en(i-1);
    end

    corBandit(i) = choice_en(i)==s(i); % correct if choice matches state
    rewdist = 1+corDots(i)+corBandit(i);
    if rewdist == 3
        count_error = 0;
    else
        count_error = count_error+1;
    end
    R(i) = round(normrnd(rewMu(rewdist),rewSigma(rewdist))); % truncate? 
        
    p = [normpdf(R(i),rewMu(1),rewSigma(1));normpdf(R(i),rewMu(2),rewSigma(2));normpdf(R(i),rewMu(3),rewSigma(3))];
    max_p = find(p==max(p));
    if max_p==3
         mue(i) = q;
         fb_p(i) = 1;
    end
    
    if ve(i-1) <= Be(end)
        ve(i) = ve(i-1)+ normrnd((mue(i)-lambda*ve(i-1)),sigma_e);
    end

    
    if count_error > 0
        h = p/sum(p)- (1-expectedPctCorr(i)); % experienced hazard rate
        h = h(max_p);
        h_s = B0 + (h^gamma/((h^gamma+(1-h)^gamma)^(1/gamma)))*(1-B0);
        h_subjective = [h_subjective h_s];
        
        if count_error >2 
            hazard_temp = hazard_temp + omega* log(1-h_subjective(count_error));
        end
        hazard = Be0 + (-log(h_subjective(1)) +log(1-h_subjective(1)) + hazard_temp);
        % not sure how to deal with the negative Be
        if hazard < 0
            hazard = 0;
        end
        Be = [Be hazard];
    else
        Be = [Be Be(end)];
        h_subjective =[] ;
        hazard_temp = 0;
       
    end
    if ve(i)<0
        ve(i)=0;
    end
end


figure
plot(ve, '*-')
%hold on
%plot(fb_p, 'o')
hold on
plot(choice_en, 'o')
hold on 
plot(s, 'k-')
hold on
%plot(1:1:ntrials, repmat(Be, ntrials, 1), '-')
plot(1:1:ntrials, Be, '-')
legend('switching variable', 'environment choice', 'environment', 'switching bound')

%% check real data

% All data is stored in 'data' structure. Rows indicate blocks, while columns indicate trial number within the blocks.
% data.choice - gives the choice made on the first of two choices ("environment/rule choice")
% -> 1 is left, 2 is right
% data.guess - gives the choice made on the second of two choices ("dots choice")
% -> 1 is left, 2 is right
% data.reward - gives the reward on any trial
% data.coherence - gives the dots stimulus coherence on any trial
% data.choicetime - gives the time (RT) it took to make the environment choice
% data.rt - gives the time it took to make the dots choice
% data.accuracy - gives the accuracy on the dots choice (right/wrong)
% -> 1 is right, 0 is wrong
% data.conf - gives the confidence reported on the dots choice
% -> can be negative or over 1, this just means they made a saccade above the bar.
% data.direction - gives the direction of the stimulus
% -> 0 is right, 180 is left
% data.highReward - gives the option that was chosen as high reward (this was once useful to me, but is not very helpful now)
% -> 1 is left, 2 is right
% data.ruleSwitch - gives the trial numbers for any one block where the environments changed
% -> the rule changes ON the trial number given
% data.stddev - gives the standard deviation of the reward distributions for a block 
% -> 1, 2 or 3. in the 'round1' dataset just deviations of 1 and 2.

clear
load('human_bandit_round1_signed.mat')

% convert to signed coh
data.coherence(data.direction==180) = -data.coherence(data.direction==180);

coh = data.coherence(:);
dotsChoice = data.guess(:)-1;
RT = data.rt(:);
conf = data.conf(:);

cohs = unique(data.coherence);
for C = 1:length(cohs)
    I = coh==cohs(C);
    pRight(C,1) = sum(I & dotsChoice==1) / sum(I);
    meanRT(C,1) = mean(RT(I));
    meanConf(C,1) = mean(conf(I));
end

figure(123); clf; set(gcf,'Position',[900   465   480   750]);
subplot(3,1,1); plot(cohs,pRight(:,1),'co-'); ylim([0 1]);
xlabel('Motion strength (%coh)'); ylabel('Proportion rightward choices');

subplot(3,1,2); plot(cohs,meanRT(:,1),'go-');
xlabel('Motion strength (%coh)'); ylabel('Reaction time (ms)');

subplot(3,1,3); plot(cohs,meanConf(:,1),'ro-'); ylim([0 1]);
xlabel('Motion strength (%coh)'); ylabel('Confidence rating');

% 
% % reconstruct reward dists
% I = data.choice==
% data.reward
    % need a vector of which bandit choice was correct



