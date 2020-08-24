%% build an experiment
clear all
close all
%mvl: Side note it be cool to check dimensionality reduction on Neurons
%group by Choice and Stim presentation (pos or negative direction) and
%compare the difference. And maybe compare the lower dimension of both
%together.
tic

% plot flags
plotMT = 0;
plotCorr = 0;

nTrials = 1000;

%seperation
seperation = 120; %30,60,90,120, and 0

% viewing duration drawn from truncated exponential distribution
dur = zeros(nTrials,1);
meandur = 360; mindur = 80; maxdur = 900;
for n = 1:nTrials
    while dur(n)<=mindur || dur(n)>=maxdur % simple trick for truncating
        dur(n) = round(exprnd(meandur));
    end
end

dur(:) = 1000; % temp: set dur to 1s (or any fixed value) to verify fano
% factor (below)

% signed motion coherence; negative is leftward
cohs = [-0.512 -0.256 -0.128 -0.064 -0.032 -1e-8 1e-8 0.032 0.064 0.128 0.256 0.512];
coh = randsample(cohs,nTrials,'true')';

dir = nan(size(coh));
dir(coh<0) = 180;
dir(coh>0) = 0;

toc

%% build MT population
tic

% nNeurons = 36;
% nNeurons = 180;
nNeurons = 360;

dirAxis = 0:360;
K = 3; % inverse variance term (will scale with coh)... mvl: Read Primer paper (might be in there, how to get K)
ampl = 60;  % actual peak FR will be about half of this, at 51.2% coh
offset = 0; % this is not the same as spontaneous firing rate, which is of 
            % course >0; keeping this offset small or zero allows driven
            % firing rates to get close to zero for nonpreferred (null)
            % direction at high coherence, as MT does

prefDirs = linspace(dirAxis(1),dirAxis(end),nNeurons+1); %mvl: Pref direction for all the neurons (One neuron for each direction...)
prefDirs(end) = [];
poscohs = cohs(cohs>0); %mvl: For positive coherences
tuning = cell(1,length(poscohs));
if plotMT
    figure(2); 
    set(gca,'xlim',[0 360],'Xtick',0:60:360);
    xlabel('direction (deg)'); 
    ylabel('firing rate (sp/s)');
    title('tuning curve for leftward preferring neuron (all cohs)');
    hold on;
end
%mvl: This is where we have to modify the tuning curve for transparent
%motion
for c = 1:length(poscohs) %mvl: Tuning curve for Each Neuron... And For Each Coherence.
    tuning{c} = zeros(nNeurons,length(dirAxis)); %mvl: Number of Neurons x Each Direction (Stimulus)... Each Neurons tuning curve. 
    for n = 1:nNeurons
        % von Mises tuning functions (captures asymmetry in pref/left response to coh)
        k = K * poscohs(c); % scale inverse variance by coh
        tuning{c}(n,:) = ampl * exp(k*cosd(dirAxis-prefDirs(n))) / (2*pi*besseli(0,k)) + offset; %mvl: Response cycle (tuning curve)
                 %^ rows are neurons, columns are dirs           
    end
    if plotMT
        plot(tuning{c}(round(nNeurons/2),:)); ylim([0 45]); hold on; 
    end
end

% set baseline to average 0% coh response (or could be a bit lower);
% this will be used below during spike train generation
baseline = mean(mean(tuning{1}));

% As an example, define two pools of neurons supplying evidence for the two
% alternatives; later can replace this with a smoothly varying weighting
% function (more realistic)
rightpool = prefDirs<60 | prefDirs>300; % within 60 deg of zero
leftpool = prefDirs>120 & prefDirs<240; % within 60 deg of 180
if plotMT
    figure(3); 
    %set(gcf,'Position', [300 300 1200 480]);
    subplot(1,2,1); 
    plot(dirAxis,tuning{6}(rightpool,:)); 
    title('tuning curves: right pool'); 
    ylim([0 45]);
    xlabel('direction (deg)'); 
    ylabel('firing rate (sp/s)');
    subplot(1,2,2); 
    plot(dirAxis,tuning{6}(leftpool,:)); 
    title('tuning curves: left pool'); 
    ylim([0 45]);
    xlabel('direction (deg)'); 
    ylabel('firing rate (sp/s)');
end

toc
%% MVL: Make a Transparent Tuning Curve with 30, 60, 90 degrees of seperation. Only for .256 Coherence. 
degOfSep = 30:30:120; %Degrees of seperation
scaling = .5; %Scaling of Seperated Direction
transTuning.S30 = cell(1,length(poscohs));
transTuning.S60 = cell(1,length(poscohs));
transTuning.S90 = cell(1,length(poscohs));
transTuning.S120 = cell(1,length(poscohs));
figure;
%construction of transparent motion tuning curves
for i = degOfSep
    for c = 1:length(tuning) 
        shiftTuning1 = [tuning{c}(:,i/2+1:end), tuning{c}(:,1:i/2)]; %Tuning curve shift for Degree of seperation in transparent stimuli, in order to scale both directions
        shiftTuning2 = [tuning{c}(:,end-i/2+1:end), tuning{c}(:,1:end-i/2)]; %Tuning curve shift for Degree of seperation in transparent stimuli, in order to scale both directions
        if i == 30
            transTuning.S30{c} = shiftTuning2.*scaling + shiftTuning1.*scaling;
            if plotMT && c == length(tuning)
                hold on;
                plot(tuning{c}(round(nNeurons/2),:), 'm--'); ylim([0 45]); 
                plot(transTuning.S30{c}(round(nNeurons/2),:), 'r'); ylim([0 45]); hold off; 
            end
        elseif i == 60 %&& c == length(tuning)
            transTuning.S60{c} = shiftTuning2.*scaling + shiftTuning1.*scaling;
            if plotMT && c == length(tuning)
                hold on;
                plot(transTuning.S60{c}(round(nNeurons/2),:), 'b'); ylim([0 45]); hold off; 
            end
        elseif i == 90 %&& c == length(tuning)
            transTuning.S90{c} = shiftTuning2.*scaling + shiftTuning1.*scaling;
            if plotMT && c == length(tuning)
                hold on;
                plot(transTuning.S90{c}(round(nNeurons/2),:), 'k'); ylim([0 45]); hold off; 
            end
        elseif i == 120 %&& c == length(tuning)
            transTuning.S120{c} = shiftTuning2.*scaling + shiftTuning1.*scaling;
            if plotMT && c == length(tuning)
                hold on;
                plot(transTuning.S120{c}(round(nNeurons/2),:), 'g'); ylim([0 45]); hold off; 
            end
        end
    end
end
title('Tuning Curve of Transparent Motion')
xlabel('Deg')
ylabel('Spike Rate')
%% make spike trains

% % independent Poisson (slow! and unrealistic)
% Counts = nan(nTrials, nNeurons);
% tic
% for t = 1:nTrials
%     c = poscohs == abs(coh(t));
%     for n = 1:nNeurons
%         lambda = tuning{c}(n,dirAxis==dir(t)) / 1000; % /1000 because tuning is in spikes/s and we're sampling every 1 ms
%         Counts(t,n) = sum(poissrnd(lambda, 1, dur(t)));
%     end
% end
% toc


% OR, spike trains according to correlation matrix


% first create corr mat
maxCorr = 0.2;
    % pairwise correlation varies smoothly from maxCorr to zero based on
    % tuning similarity (pref dir proximity)
fano = 1.8; % fano factor = variance/mean ratio (forced to be 1 for poisson,
            % but for real neurons it's often 1.5-1.8)
c=linspace(maxCorr,0.0,nNeurons/2); 
c=[c fliplr(c)];
Cor=toeplitz(c); % aka diagonal-constant matrix, shortcut for making the kind of corr mat we want
Cor(logical(eye(size(Cor)))) = 1; % set main diagonal to be 1
tmpCor = Cor; 
tmpCor(tmpCor==1)=NaN; % for better color range when plotting
if plotCorr
    figure(4); imagesc(tmpCor); set(gca,'ydir','normal'); axis square; colorbar;
    title('spike count correlation matrix: intended');
    xlabel('neuron ID (pref dir)'); ylabel('neuron ID (pref dir)');
end


% then use corr mat to generate spike  counts across the population.
% here are two (old) methods; there may be better ones more recently

% % (1) Mazurek et al. 2002 Nat Neurosci (intuitive, but slower)
% % % % for g = 1:20 % TEMP: to test that average corr value is as desired
% Counts = nan(nTrials, nNeurons);
% tic
% for t = 1:nTrials
%     i = dir(t);
%     c = poscohs == abs(coh(t));
%     m = tuning{c}(:,dirAxis==dir(t)); % mean
%     T = dur(t)/1000;
%         % now we have a vector of mean counts and a correlation matrix based on
%         % tuning similarity. So to get a covariance matrix just multiply by the
%         % product of ith and jth s.d. where s.d. is sqrt(fano*m*T)
%     sd = sqrt(fano*m*T);
%     Cov = repmat(sd,1,nNeurons) .* repmat(sd,1,nNeurons)' .* Cor;
%     Counts(t,:) = round(mvnrnd(m*T, Cov)); % multivariate normal random numbers
% end
% Counts(Counts<0)=0;
% toc
% % % % % plot pairwise responses for similar and unsimilar neurons, for a fixed stimulus (or not)
% % % % whichTr = abs(coh)<0.001;
% % % % durToNorm = dur(whichTr);
% % % % corrs1(g) = corr(Counts(whichTr,1)./durToNorm, Counts(whichTr,2)./durToNorm);
% % % % corrs2(g) = corr(Counts(whichTr,1)./durToNorm, Counts(whichTr,19)./durToNorm);
% % % % % figure; plot(Counts(whichTr,1),Counts(whichTr,2),'x'); title(num2str(corrs1(g)));
% % % % % figure; plot(Counts(whichTr,1),Counts(whichTr,19),'o'); title(num2str(corrs2(g)));
% % % % end
% % % % mean(corrs1)
% % % % mean(corrs2) % TEMP: to test that average corr value is as desired



% (2) Shadlen et al. 1996 J Neurosci, Appendix 1 (>2x faster)
% Change with tuning curve to use
%decide which tuning curve we using, and if using transparent motion choose
%seperation
if seperation == 0
    tuningcurve = tuning;
elseif seperation == 30
    tuningcurve = transTuning.S30;
elseif seperation == 60
    tuningcurve = transTuning.S60;
elseif seperation == 90
    tuningcurve = transTuning.S90;
elseif seperation == 120
    tuningcurve = transTuning.S120;
end
tic
rootQ = sqrtm(Cor);
Counts = nan(nTrials, nNeurons);
for t = 1:nTrials
    c = poscohs==abs(coh(t)); %mvl: Which coherence are we working with (regardless of sign) in this trial
    T = dur(t)/1000; %mvl: Turn into second(s) instead of millisecond(ms)
    m = tuningcurve{c}(:,dirAxis==dir(t)); %tuning{c}(:,dirAxis==dir(t)); 
    %mvl: Which direction are we working with (including coherence), use that to cluster the mean rate of all the neurons for particular direction
    %Population Curves will be when TransMotion is motion evenly between
    %180 or 0 degrees.
    z = normrnd(0,1,nNeurons,1); %mvl: Add noise to each neuron (360x1)
    y = rootQ * z; %mvl: The correlation matrix acts as weights to attribute the correct proportion of effect X neurons have on the Y neuron, and that makes up the scaliability that is them multiply to the variance 
    y = y.*sqrt(fano*m*T) + m*T; % variance is fano*mean/second, and SD is sqrt of that %mvl: Response with Noise (Noise = CorrNoise*IndNoise), mean and variance are measurements of time (per seconds)
    %MVL: basicalyl the first part of the top equation (before the '+') is
    %normrnd(0, sqrtm(Cor).*sqrt(fano*m*T), nNeurons, 1). Each neuron has a
    %variance (fano*m*T). 
    y(y<0) = 0; % no negative reponses
    Counts(t,:) = round(y); %mvl: Trial x Neurons Firing Counts (A noisy rendition, instead of the mean vonmiss)
end
toc


% below we'll also need spike *rates*, since dur varies
Rates = Counts./(dur/1000); %mvl: Convert Firing 'counts' to rates (spikes/second)

%% check how well we achieved the intended corr mat & fano

if plotCorr 
    figure(5); 
    %set(gcf,'position',[520 200 1500 1100]); 
end
Fanos = nan(nNeurons,length(cohs));
for c = 1:length(cohs)
    I = coh==cohs(c);
    ratemat = Rates(I,:); %mvl: Select FiringRates [Spikes/Second] (based on trials) that belong to this Coherence (All neurons in those trials)
    corrmat = corrcov(cov(ratemat));
    
%     % uncomment this to reassure yourself that corrcov() works:
%     if c==1
%     corrmat2 = nan(size(corrmat));
%     for i = 1:nNeurons
%         for j = 1:nNeurons
%             corrmat2(i,j) = corr(ratemat(:,i),ratemat(:,j));
%         end
%     end
%     figure;plot(corrmat(:),corrmat2(:),'x');
%     end

    if plotCorr
        figure(5);
        subplot(3,4,c); 
        imagesc(corrmat); 
        set(gca,'ydir','normal'); 
        axis square; 
        colorbar;
        xlabel('neuron ID (pref dir)'); 
        ylabel('neuron ID (pref dir)');
        caxis([-0.1 maxCorr+0.1]); title(num2str(cohs(c)));
    end
    
    % variance to mean ratio across repeated presentations of the same
    % stimulus for all neurons; will be averaged across stimuli (cohs)
    % mvl: Its not the same stimulus just the same Coherence (I guess we
    % can say this)
    Fanos(:,c) = var(Counts(I,:))./mean(Counts(I,:)); %mvl: Fano of each Neuron 
%     Fanos(:,c) = var(Rates(I,:))./mean(Rates(I,:));
% **this greatly overestimates the true Fano factor, even when using Rates, 
% because of added variance from variable durations. set dur(:)=1000 on 
% line ~23 to see the correct Fano

end
if plotCorr
    mean(Fanos)
    mean(mean(Fanos))
    figure(5); %suptitle('spike count correlation matrix: actual (sep by coh)');
end

    
% confirm that variance of difference variable (momentary evidence)
% increases with coh, aka signal-dependent noise
if plotMT
rDiff = nan(nTrials,1); %mvl: Difference in Firing Rates between population
for t = 1:nTrials
    rDiff(t) = sum(Counts(t,rightpool)) - sum(Counts(t,leftpool));
end
for c = 1:length(poscohs)
    I = abs(coh)==poscohs(c);
    rDiffMean(c) = mean(abs(rDiff(I))); %mvl: Mean in difference in Population (Right vs Left)
    rDiffVar(c) = var(rDiff(I)); %mvl: Variance in differene in Population (Right vs Left)
end
figure(6); 
subplot(1,2,1); 
%set(gcf,'position',[500 500 1000 420]);
plot(poscohs,rDiffMean,'o-'); 
title('mean of difference variable (momentary evidence)'); 
xlabel('unsigned coh');
subplot(1,2,2); 
plot(poscohs,rDiffVar,'v-'); 
title('variance of difference variable (momentary evidence)'); 
xlabel('unsigned coh');
end


%% Back to Mazurek for spike train generation (inverse probability method):

% construct spike density function for each trial according to known time
% course of MT responses, convert to cumulative probability (normalized),
% and use to map uniform random deviates (one for each spike in Counts)
% onto spike times for every trial

% an interesting exercise (which I've yet to do) is check whether this
% method generates the expected cross-correlograms (see e.g. Bair et
% al. J Neurosci 2001). Could be that the 'dummy spike train' from 
% Mazurek is needed for this level of realism.


% baseline = 10; % can be specified independently of tuning curves, but should
               % be higher than response to high-coh motion in antipreferred
               % (null) direction;
       % for now this is set above based on tuning{}               
latency = 50; % ms
risetime = 20; % ms
dT = 1;
tAxis = 1:dT:max(dur)+3*latency;
               % it's 3x because 1x at the beginning, 2x at the end (ramp-down is slower than rise time)

R = nan(nTrials, nNeurons, length(tAxis)); %mvl: 3D Mat with Time being 3rd Dimension
disp('Generating spike trains, may take tens of seconds to minutes...');
tic
for t = 1:nTrials
    for n = 1:nNeurons %mvl: iterated through each trial and then neuron
                
        meanrate = Rates(t,n);

        sdf=[]; % spike density function
        
            % lay down some baseline rate, from t=1 to vis latency
        sdf(1:latency) = baseline;
            % quick linear rise to mean
        sdf(end+1 : latency+risetime) = linspace(baseline,meanrate,risetime);
%             % steady state with some decay*, like real MT
            decay = 1+0.2*sign(baseline-meanrate); % *up or down depending on whether mean rate is above or below baseline
        sdf(end+1 : latency+dur(t)) = linspace(meanrate, meanrate*decay, dur(t)-risetime);
              % include die-down over 2x latency
        sdf(end+1 : end+2*latency) = linspace(meanrate*decay, baseline, 2*latency);

% %         % check what it looks like, if you want
% %         keyboard
% %         figure; plot(sdf); % CAREFUL! don't run entire loop when uncommented!
                
        % The Mazurek method (I think) assumes zero noise in the driving
        % input (underlying 'rate' or mean of poisson process, see
        % Churchland et al. 2011 Neuron paper "Variance 1as a Signature...")
        % which is unrealistic even for MT. At the very least it is
        % a function of motion energy, plus any noise added between retina
        % and MT. So adding some ad hoc variability is easily justified
        sdfNoiseSD = 0.2 * sdf; % arbitrary, but scales with spike rate
            % make fluctuations on the time scale of frames (motion energy)?
        % SDFnoise(1:8:dur(t)+2*latency) = SDFnoiseSD*randn(1,length(1:8:dur(t)+2*latency));
            % no, for now just every ms; could use actual trials' ME later
        sdf = sdf + sdfNoiseSD.*randn(1,length(sdf));
        sdf(sdf<0)=0;
        
% %         % may also wish to inspect the noisy version
% %         keyboard
% %         figure; plot(sdf); % CAREFUL! don't run entire loop when uncommented!
                            
        csdf = cumsum(sdf); % cumulative spike density
        cspf = csdf / csdf(end); % cumulative spike probability
            % now we simply draw nspikes uniform random numbers on [0 1] 
            % and use the CSPF as a lookup table for where to place them in
            % time*, redrawing any that fall within a refractory period
            % [*see note below]

            % interpolate to ensure intersect
            % cspfI = interp1(1:length(cspf),cspf,1:0.1:length(cspf));  ??? don't remember what this was about
        
        R(t,n,1:length(sdf)) = 0; % initialize spike train
        % ^ trial, neuron, time    

        sptimes = nan(1,Counts(t,n));
        for s = 1:Counts(t,n) % is there a faster way than spike by spike?
            rnd = rand;
            sp = find(abs(cspf-rnd)==min(abs(cspf-rnd)));
            sp = sp(1); % avoid non-unique times
            if s>1
                while any(abs(sptimes-sp)<2) %mvl: Dont allow spikes to go off 2ms after one another, redraw is so
                % ^ enforce refractory period by redrawing rnd
                    rnd = rand;
                    sp = find(abs(cspf-rnd)==min(abs(cspf-rnd)));
                    sp = sp(1); % avoid non-unique times
                end
            end
            sptimes(s) = sp;
        end
        R(t,n,sptimes) = 1;
        
    end
end
disp('done.');
toc


% check that we achieved the intended spike counts
if plotMT
    Counts2 = nansum(R,3);
    figure(7); 
    plot(Counts(:),Counts2(:),'x');
    title('intended vs actual spike counts (should be identical)');
end

% * Note that this method distributes the allocated spikes (determined by
% tuning curves and correlations) throughout the trial, which now
% includes a baseline and ramp-down period. In real experiments, tuning
% curves are usually constructed by counting spikes in a window that
% excludes these periods. This has two consequences: (1) spike counts/rates
% in the appropriate time windows of these constructed spike trains will be
% lower than the intended values based on tuning{}, and (2) the average
% spike rate during the baseline period (see PSTHs below) will depend on
% the mean rate for that trial, i.e. will depend on coherence, which of
% course would not happen in a real experiment. I don't think this should
% cause any problems with this exercise but noting it here just in case.



%% some sanity checks:

if plotMT

% PSTHs, do they look reasonable?
cohsForPSTH = [-0.512 -0.128 -0.032 0.032 0.128 0.512];
clr = cool(length(cohsForPSTH));
figure(8); 
for c = 1:length(cohsForPSTH)
    I = coh==cohsForPSTH(c); % R = [trials, neurons, time]
    plot(nanmean(squeeze(nanmean(R(I,rightpool,1:800),2)))*1000,'Color', clr(c,:)); hold on;
end
title('right pool'); 
xlabel('time from dots onset (ms)'); 
ylabel('spikes/s');
legend(cellstr(num2str(cohsForPSTH')),'location','north','orientation','horizontal');

figure(9); 
for c = 1:length(cohsForPSTH)
    I = coh==cohsForPSTH(c);          % R = [trials, neurons, time]
    plot(nanmean(squeeze(nanmean(R(I,leftpool,1:800),2)))*1000,'Color', clr(c,:)); hold on;    
end
title('left pool'); 
xlabel('time from dots onset (ms)'); 
ylabel('spikes/s');
legend(cellstr(num2str(cohsForPSTH')),'location','north','orientation','horizontal');


% what about raster plots?
figure(10);
C = coh==0.512;
raster = squeeze(R(C,1,1:400));
[I,J] = find(raster>0);
plot(J,I,'.'); 
xlabel('time from dots onset (ms)'); 
ylabel('trial'); 
ylim([0 size(raster,1)]);
title('spike raster: across trials for a given neuron');
 
figure(11);
C = find(coh==0.512);
raster = squeeze(R(C(1),:,1:400));
[I,J] = find(raster>0);
plot(J,I,'.'); 
xlabel('time from dots onset (ms)'); 
ylabel('neuron (pref dir)'); ylim([0 nNeurons]);
title('spike raster: across neurons for a given trial');

end
        

%% save the results
disp('saving...');

if seperation == 0
    clearvars -except R tuning nNeurons nTrials prefDirs latency dir dur coh cohs poscohs tAxis
    tic
    save(sprintf('simMT_nNeu=%d_nTr=%d.mat', nNeurons, nTrials),'-v7.3');
    toc
else %seperation leads to Transparent Tuning
    clearvars -except R tuning nNeurons nTrials prefDirs latency dir dur coh cohs poscohs tAxis transTuning seperation
    tic
    save(sprintf('TransMotion=%d simMT_nNeu=%d_nTr=%d.mat', seperation, nNeurons, nTrials),'-v7.3');
    toc
end

disp('done.');




