% best 6 param (dsigma2 only) = 35588.44 (run3 in home dir)
% k= 0.283463
% B= 31.971901
% theta= 0.603722
% dcoh?stim= 0.131371
% dcoh?bias= ?0.013896
% dk= 0.000000
% dsigma= 0.305417
% sigma2?beta= 0.000000
% weber= 0.000000
% guessrate= 0.000000
% theta2= theta
% dtheta= 0.000000
% dtheta2= 0.000000


% best 6 param (yoked dthetas only) = 35648.266962 (run14) [confirmed 8x]

% best 7 param (dsigma2+yoked dthetas) = 35574.909128 (screenshot)
    % BUT allows dsigma and dtheta to work in opposite directions, leading to large
    % values of dsigma -- not very parsimonious
    
% best 7 param w/ dtheta constrained negative: 35591.272197 - CONFIRMED GLOBAL
% ****************************************
% run 2
% 	k= 0.281146
% 	B= 31.977
% 	theta= 0.596453
% 	dcohstim= 0.130977
% 	dcohbias= -0.0137389
% 	dk= 0
% 	dsigma2= 0.274418
% 	sigma2beta= 0
% 	weber= 0.000000
% 	guessrate= 0
% 	theta2= 0.596453
% 	dtheta= -0.0014188
% 	dtheta2= -0.0014188
% err: 35591.272197
% otherwise, if dtheta is started higher and dsigma lower, always hits local min, e.g. run 5 on home dir (35615)


% best 7 param (dsigma2 and SIGMA2BETA FTW!) = 35556.536347 (run13)  ***

% best 7 param (sigma2beta and dthetas) = 35605.111491 (run 15) [confirmed 4x]

% 
% LLnull = -35588.44;
% df = 6;
% n = 53134;
% BICnull = -2*LLnull + df*log(n)
% 
% LLalt = -35648.266962;
% df = 7;
% n = 53134;
% BICalt = -2*LLalt + df*log(n)
% 
% deltaBIC = BICalt - BICnull
% AltSigBetter = (BICalt < BICnull)



clear all
close all

% stim_type = 'pseudostim';
% stim_type = 'elestim';
% stim_type = 'optstim';
stim_type = 'optstimB';

% filename = 'FIRA_uStim_bothmonks_All_N=63'; PsMag = NaN;
% filename = 'FIRA_uStim_Ike_All_N=39'; PsMag = NaN;
% filename = 'FIRA_uStim_Damien_All_N=24'; PsMag = NaN;
% filename = 'FIRA_PsStim_All_Sessions_N=42'; PsMag = NaN; stim_type = 'pseudostim'; % g_stimeff = PsMag/100;
% filename = 'FIRA_uStim_Damien_highcurr_N=5'; PsMag = NaN;
% filename = 'workspace_NoSlopeReduction_N=48'; PsMag = NaN;
% filename = 'FIRA_uStim_Damien_highcurr_N=8'; PsMag = NaN;

folder = '/Users/crfetsch/Documents/MATLAB/projects/Opto/';
% filename = 'FIRA_jaws';
filename = 'jaws_subset_forFitting'; 
PsMag = NaN;



savePDF = 1;
saveMAT = 0;
fourCurves = 0;
paper_figs = 1; %CAREFUL: MAY REPLACE FIGS IN ILLUSTRATOR FILE

% load([folder filename '.mat'], 'FIRA');
    % TEMP: use the all_ind leftover from jaws datamining
load([folder filename '.mat'], 'FIRA', 'all_ind');

stimtrg=2;

    %extract some of the relevant parameters from FIRA
% all_ind = getTrialIndex_certstim(FIRA,[],[],[],[0 inf],[],[0 1 2]);    %all trials
    % TEMP: use the all_ind leftover from jaws datamining
coh_ = FIRA{2}(all_ind,IND('dot_coh'));
dur_ = FIRA{2}(all_ind,IND('dot_off')) - FIRA{2}(all_ind,IND('dots_on'));
delay_ = FIRA{2}(all_ind,IND('fp_off')) - FIRA{2}(all_ind,IND('dot_off'));
strg_ = ~isnan(FIRA{2}(all_ind,getFIRA_columnByName('s_trg')));
cor_ = FIRA{2}(all_ind,IND('correct'));
cor_trg_ = FIRA{2}(all_ind,IND('cor_trg'));
cho_trg_ = FIRA{2}(all_ind,IND('cho_trg'));
coh_(cor_trg_~=stimtrg) = -coh_(cor_trg_~=stimtrg);
stim_ = FIRA{2}(all_ind,IND(stim_type));
stim_(isnan(stim_)) = 0;
% session_ = FIRA{1}.ownership;
coh_set = {-512, -256, -128, -64, [-32 -16], 0, [16 32], 64, 128, 256, 512};  %how to group coherences
coh_axis = [-512, -256, -128, -64, -32, 0, 32, 64, 128, 256, 512];      %what coherence is representative of each group
coh_set_unsigned = coh_set(6:end);
coh_axis_unsigned = coh_axis(6:end);
dur_(dur_<60) = 60; % kluge a strange zero-dur trial

    
    %fit Pright and Psure of non-stimulated trials and Pright of stimulated trials to calculate k, B,
    %logodds_crit and init_bias, establishing a *prediction* for Psure of stim trials
% FitStrategy = 'OneStage-PrightPS,Pright'; strategyIndex = 'A';

    % fit Pright and Psure of all trials to calculate k, B, logodds_crit and init_bias
FitStrategy = 'OneStage-PrightPS,PrightPS'; strategyIndex = 'B';

% FitStrategy = 'TwoStage-PrightPS,PrightPS'; strategyIndex = 'B';


    %run the fitting procedure for the specified strategy
switch FitStrategy,
    case {'TwoStage-PrightPS,Pright','TwoStage-PrightPS,PrightPS'},
        if strcmp(FitStrategy,'TwoStage-PrightPS,Pright'),
            fitwhat = {'PrightPS', 'Pright'}; 
        else
            fitwhat = {'PrightPS', 'PrightPS'}; 
        end;
            %first stage: fit Pright and Psure of non-stimulated trials to calculate k, B, theta and init_bias
%        1  2 3     4      5    6  7     8        9      10        11     12     13
%       [k  B theta stmeff bias dk dsig2 sig2beta weber guess_rate theta2 dtheta dtheta2] %%%%%%
fixed = [0  0 0     0      0    1  0     1        1      1         1      1      1];
guess = [
    0.294
    31.6
    0.609
    0.112
    -0.0134
    0
    0.237
    0
    0
    3.7201e-44
    0.609
    0
    0]';

        L = find(stim_==0);
        D.cor_trg    = cor_trg_(L);
        D.cho_trg    = cho_trg_(L);
        D.coh     = coh_(L)/1000;      %must be signed coherence, this is different from the Science paper 
        D.dur = dur_(L);
        D.SureT   = strg_(L);
        D.ustim = stim_(L);
            clear options;
            options.fitwhat = fitwhat;
            options.yoke_theta = 1;
            options.fitMethod = 'fms';
        [modelParam(1),modelLL(1),trialData(1)] = fit_Diffusion_posterior_new4(D, guess, fixed, 1, options);
        expectedPright{1} = trialData(1).expectedPright;
        expectedPS{1} = trialData(1).expectedPS;

    %second stage: fit stimulated trials to calculate uStim_eff and thetaOffset
%        1  2 3     4      5    6  7    8        9      10        11     12     13
%       [k  B theta stmeff bias dk dsig sig2beta weber guess_rate theta2 dtheta dtheta2] %%%%%%
fixed = [1  1 1     0      1    0  0    0        1      1         0      0      0];
guess = [
    modelParam.final(1)
    modelParam.final(2)
    modelParam.final(3)
    0.124
    modelParam.final(5)
    -0.05
    0.22
    0
    0
    3.7201e-44
    0.578441
    0
    0]';

        L = find(stim_==1);
        D.cor_trg    = cor_trg_(L);
        D.cho_trg    = cho_trg_(L);
        D.coh     = coh_(L)/1000;      %must be signed coherence, this is different from the Science paper 
        D.dur = dur_(L);
        D.SureT   = strg_(L);
        D.ustim = stim_(L);
            clear options;
            options.fitwhat = fitwhat;
            options.yoke_theta = 1;   
            options.fitMethod = 'fms';
        [modelParam(2),modelLL(2),trialData(2)] = fit_Diffusion_posterior_new4(D, guess, fixed, 1, options);
        expectedPright{2} = trialData(2).expectedPright;
        expectedPS{2} = trialData(2).expectedPS;
    case {'OneStage-PrightPS,Pright','OneStage-PrightPS,PrightPS','OneStage-PS,PS','OneStage-PrightPS,Pright-woTs'},
        if strcmp(FitStrategy,'OneStage-PrightPS,Pright'),
            %fit Pright and Psure of non-stimulated trials and Pright of stimulated
            %trials to calculate k, B, logodds_crit and init_bias, 
            %establishing a PREDICTION for Psure of stimulated trials
            fitwhat = {'PrightPS', 'Pright'};
        elseif strcmp(FitStrategy,'OneStage-PrightPS,PrightPS'),
            %fit Pright and Psure of all trials to calculate k, B, logodds_crit and init_bias
            fitwhat = {'PrightPS', 'PrightPS'}; 
        elseif strcmp(FitStrategy,'OneStage-PS,PS'),
            fitwhat = {'PS', 'PS'};
        else
            fitwhat = {'PrightPS', 'Pright-woTs'};
        end;

% %         load fitParam_6-27-13_SigFree_CohVar_FIT;
% %         load fitParam_6-30-13_SigFree_CohVar+km_FIT; guess=fitParam; % 333772.9050
%         load fitParam_6-30-13_SigFree_CohVar+km_PRED; guess=fitParam;

%         guess(2)=30;
%         guess(4)=[];
%         guess(6)=guess(1)*guess(8)-guess(1);
%         guess(7)=sqrt(guess(9)+1)-1;
%         guess(8)=guess(11);  %0.441264;
%         guess(9)=0;
%         guess(10)=3.7201e-44;
%         guess(11)=guess(3);
%         guess(12:13)=0;

%        1  2 3     4      5    6  7    8        9      10        11     12     13
%       [k  B theta stmeff bias dk dsig2 sig2beta weber guess_rate theta2 dtheta dtheta2] %%%%%%
fixed = [0  0 0     0      0    1  1     1        1     1          1      1      1];

% k= 0.283463
% B= 31.971901
% theta= 0.603722
% dcohstim= 0.131371
% dcohbias= -0.013896
% dk= 0.000000
% dsigma2= 0.305417
% sigma2beta= 0.000000
% weber= 0.000000
% guessrate= 0.000000
% theta2= theta
% dtheta= 0.000000
% dtheta2= 0.000000

    
% jaws
k= 0.357002
B= 22.3505
theta= 1.5 % 1.47555
dcohstim= -0.051 % -0.0378313
dcohbias= 0.0102631
dk= 0
dsigma2= 0
sigma2beta= 0
weber= 0.000000
guessrate= 0
theta2= theta
dtheta= 0
dtheta2= 0


guess = [k B theta dcohstim dcohbias dk dsigma2 sigma2beta weber guessrate theta2 dtheta dtheta2];

    fixed(:)=1; % TEMP: for hand-tuning, to make better guesses, etc
    

% find best LL for pseudostim
% w/o sig2beta:
% (table 1)

% w sig2beta:
% run1 - 28155.091882 (dsigma2 started and ended up near zero)
% run2 - 28150.641852 (dsigma2 near 0.05)
% run3 - 28154.436637 (only dsigma2 and beta free, former gets near 0.09, latter near zero)

% % best: (confirmed w/ multi)
% ****************************************
% run 2
% 	k= 0.312476
% 	B= 31.8858
% 	theta= 0.497495
% 	dcohstim= 0.164878
% 	dcohbias= -0.0229172
% 	dk= 0
% 	dsigma2= 0.096382
% 	sigma2beta= 0.427636
% 	weber= 0.000000
% 	guessrate= 0
% 	theta2= 0.497495
% 	dtheta= 0
% 	dtheta2= 0
% err: 28150.322512


% load finalFits_uStim fitParam; % 35682 w/ wrong dsigma2
%     fitParam(7) = 0.237; % 35635   w/ correct dsigma2
% guess=fitParam;
    
% load guess_LowCurSimple_Jan24.mat; % RK and I both get 35601.89 [not as good as run3 mentioned above!]
% load guess_LowCurExtended_Jan24.mat; % should be 35540.22 but I get 35542.22
% load guess_LowCurHyperextended_Jan24.mat; % should be 35541.41 but I get 35542.83
   
%     guess = guess .* (rand(size(guess))/2.5-0.2+1);
%     guess(11) = guess(3);
%     guess(13) = guess(12);
        
        allfixed = all(fixed==1);
 
        D.cor_trg = cor_trg_;
        D.cho_trg = cho_trg_;
        D.coh     = coh_/1000;      %must be signed coherence, this is different from the Science paper 
        D.dur = dur_;
        D.SureT   = strg_;
        D.ustim = stim_;
        
            clear options;
            options.fitwhat = fitwhat;
            options.yoke_theta = 0; % yoke thetas and dthetas (formerly "flag")
            options.fitMethod = 'fms';
%             options.fitMethod = 'global';
%             options.fitMethod = 'multi';
%             options.fitMethod = 'pattern';
        [modelParam(2),modelLL(2),trialData(2)] = fit_Diffusion_posterior_new4(D, guess, fixed, 0, options);
        expectedPright{2} = trialData(1).expectedPright;
        expectedPS{2} = trialData(1).expectedPS;
end;

%********** End, fit Pcor and Psure *********


%********** Begin, Make smooth Pright and Psure curves based on model parameters *********

if allfixed
    fitParam = guess;
    modelParam(2).final = fitParam;
else
    fitParam = modelParam(2).final;
end
global err2



% temp: hand-tune params -- comment out this line (cell begin) to "HERE RESUMES NORMAL CODE", below, to resume normal operation
% codebank(47) -- now can be handled with "allfixed"

% save finalFits_temp fitParam dur_ FIRA trialData

    %generate a set of random stimulus durations.  
dursamples = 2000;      % the actual number of simulated trials will be
                        % dursamples * length(g_coh) * 4
                        % the 4 comes from 2[for_Ts] * 2[for_uStim]

if ~exist('dur_','var')
    dur_ = trialData(2).dur;
end
                        
% Q(1) = 0; Q(5) = inf; Q(2:4) = quantile(dur_,0.25:0.2:0.75); % quartiles
% dur_setQ{1} = [Q(1) Q(2)];
% dur_setQ{2} = [Q(2) Q(3)];
% dur_setQ{3} = [Q(3) Q(4)];
% dur_setQ{4} = [Q(4) Q(5)];
% dur_setQ{5} = [0 inf];

% OR: median split
Q(1) = 0; Q(3) = inf; Q(2) = median(dur_);
dur_setQ{1} = [Q(1) Q(2)];
dur_setQ{2} = [Q(2)+1 Q(3)];
dur_setQ{3} = [0 inf];

% for q = 1:5 % run/plot for 4 dur quantiles as well as full dataset
for q = 3 % run/plot for 2 halves as well as full dataset (q=3 means all trials)

if q==3 %5
    rand_dur = randsample(dur_, dursamples, 'true');
else
    rand_dur = randsample(dur_(dur_>Q(q) & dur_<Q(q+1)), dursamples, 'true');
end

g_coh = unique([(-0.6:0.04:0.6)'; cat(2,coh_set{:})'/1000]);

D.coh   = repmat(g_coh', [length(rand_dur)*4 1]);
% D.dur   = repmat(repmat(rand_dur',[1 length(g_coh)]), [4 1]);  WTF?
D.dur   = repmat(rand_dur,4,length(g_coh));
D.cor_trg      = ones(size(D.coh));
D.cor_trg(D.coh<0) = 2;
D.cho_trg   = D.cor_trg;
D.SureT     = ones(size(D.coh));
D.SureT(1:2*dursamples,:) = 0;
D.ustim = ones(size(D.coh));
D.ustim([1:dursamples, 2*dursamples+(1:dursamples)],:) = 0;
    %make them vertical vectors
D.coh = D.coh(:);
D.dur = D.dur(:);
D.cor_trg = D.cor_trg(:);
D.cho_trg = D.cho_trg(:);
D.SureT = D.SureT(:);
D.ustim = D.ustim(:);

options.coh_set = [-0.512 -0.256 -0.128 -0.064 -0.032 0 0.032 0.064 0.128 0.256 0.512]';
options.coh_set_freq = [1/12*ones(1,5), 1/12*2, 1/12*ones(1,5)]';
options.coh_set_freq = [options.coh_set_freq options.coh_set_freq];

% % this is weird; why would I have set coh_set_freq to the frequencies from
% % the real data, when this is just plotting the fake data-based fit curves???
% options.coh_set = unique(trialData(2).coh);
% options.coh_set_freq = nan(length(coh_set), 2);        
% for s = 0 : 1
%     for c = 1 : length(coh_set)
%         options.coh_set_freq(c,s+1) = sum(trialData(2).coh==options.coh_set(c)&trialData(2).ustim==s)/length(trialData(2).coh);
%     end
% end

options.plot = 0;
% options.yoke_theta = 1;

[m1,m2,D] = fit_Diffusion_posterior_new4(D, fitParam, ones(size(fitParam)), 0, options);

g_pright_nostrg = nan(length(g_coh), 2);
g_pright_strg = nan(length(g_coh), 2);
g_pright_all = nan(length(g_coh), 2);
g_pstrg = nan(length(g_coh), 2);
for s = 0 : 1
    for c = 1 : length(g_coh)
        I = D.coh==g_coh(c) & D.ustim==s & D.SureT==0;
        g_pright_nostrg(c,s+1) = nanmean(D.expectedPright(I));
        I = D.coh==g_coh(c) & D.ustim==s & D.SureT==1;
        g_pright_strg(c,s+1) = nanmean(D.expectedPright(I));
        g_pstrg(c,s+1) = nanmean(D.expectedPS(I));
        I = D.coh==g_coh(c) & D.ustim==s;
        g_pright_all(c,s+1) = nanmean(D.expectedPright(I));
    end
end

%********** End, Make smooth Pright and Psure curves based on model parameters *********

    
    %********** Begin, measure Pcor and Psure *********
    all_ind = getTrialIndex_certstim(FIRA,[],[],[],dur_setQ{q},[],[0 1 2]); 
    stimtrg = 2;
    coh_set = {-512, -256, -128, -64, [-32 -16], 0, [16 32], 64, 128, 256, 512};  %how to group coherences
    coh_axis = [-512, -256, -128, -64, -32, 0, 32, 64, 128, 256, 512];      %what coherence is representative of each group
    coh_ = FIRA{2}(all_ind,IND('dot_coh'));
    dur_ = FIRA{2}(all_ind,IND('dot_off')) - FIRA{2}(all_ind,IND('dots_on'));       
%     strg_delay_ = FIRA{2}(all_ind,IND('s_on')) - FIRA{2}(all_ind,IND('dot_off'));
        strg_ = ~isnan(FIRA{2}(all_ind,getFIRA_columnByName('s_trg')));
    cor_ = FIRA{2}(all_ind,IND('correct'));                                         
    cor_trg_ = FIRA{2}(all_ind,IND('cor_trg'));
    cho_trg_ = FIRA{2}(all_ind,IND('cho_trg'));
    coh_(cor_trg_~=stimtrg) = -coh_(cor_trg_~=stimtrg);
    if ~isempty(FIRA{2}(all_ind,IND('pseudostim'))) && sum(FIRA{2}(all_ind,IND('pseudostim'))) > 5
        stim_ = FIRA{2}(all_ind,IND('pseudostim'));
    else
        stim_ = FIRA{2}(all_ind,IND(stim_type));
    end
    stim_(isnan(stim_)) = 0;

        %now explor the monkey's performance, calculate p(right), p(strg), etc
    pright_strg = nan(length(coh_set), 2);
    pcorr_strg = nan(length(coh_set), 2);
    n_strg = nan(size(pright_strg));
    pright_nostrg = nan(length(coh_set), 2);
    pcorr_nostrg = nan(length(coh_set), 2);
    n_nostrg = nan(size(pright_nostrg));
    pstrg = nan(length(coh_set), 2);
    nstrg = nan(size(pstrg));
    pright_all = nan(length(coh_set), 2);


    % J = ~isnan(strg_delay_);    %find trials with sure target
    % bhLog has no strg_delay. use suitable replacement:
    J = ~isnan(FIRA{2}(all_ind,getFIRA_columnByName('s_trg')));    %find trials with sure target

    %calculate Pright and Pstrg for signed coherences 
for s = 0 : 1,
    S = stim_==s;
    I = J & S & ismember(cho_trg_,[1 2]);
    [pright_strg(:,s+1), pright_strg_se(:,s+1)] = calcGroupMean(cho_trg_(I)==stimtrg, coh_(I), coh_set, 'binary');
    I = J & S;
    [pstrg(:,s+1), pstrg_se(:,s+1)] = calcGroupMean(cho_trg_(I)==3, coh_(I), coh_set, 'binary');
    I = ~J & S & ismember(cho_trg_,[1 2]);
    [pright_nostrg(:,s+1), pright_nostrg_se(:,s+1)] = calcGroupMean(cho_trg_(I)==stimtrg, coh_(I), coh_set, 'binary');
    I = S & ismember(cho_trg_,[1 2]);
    [pright_all(:,s+1), pright_all_se(:,s+1)] = calcGroupMean(cho_trg_(I)==stimtrg, coh_(I), coh_set, 'binary');
end;


%********** Begin, explore how well the fits corresponds with the monkey's behavior *********

%    Pright and Psure as a function of signed coherence 
figure(111+q); clf;
set(gcf, 'Color', [1 1 1], 'Position', [300+100*q 300 600 850], 'PaperPositionMode', 'auto');
        %Psure
subplot(3,1,1);
hold on;
errorbar(coh_axis/10, pstrg(:,1), pstrg_se(:,1), 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
errorbar(coh_axis/10, pstrg(:,2), pstrg_se(:,2), 'd', 'Color', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 4);
plot(g_coh*100, g_pstrg(:,1), 'k-');
    plot([g_coh(g_pstrg(:,1)==max(g_pstrg(:,1)))*100 g_coh(g_pstrg(:,1)==max(g_pstrg(:,1)))*100], [0 max(g_pstrg(:,1))], 'g-', 'LineWidth', 2);
plot(g_coh*100, g_pstrg(:,2), 'k-.');
    plot([g_coh(g_pstrg(:,2)==max(g_pstrg(:,2)))*100 g_coh(g_pstrg(:,2)==max(g_pstrg(:,2)))*100], [0 max(g_pstrg(:,2))], 'g-', 'LineWidth', 2);
    plot([g_coh(1) g_coh(end)]*100,[0.5 0.5], 'k--');
plot([0 0],[0 1],'k--');
set(gca, 'XLim', [-60 60], 'XTick', -60:10:60, 'XTickLabel', makeTickLabel(-60:10:60,20), ...
         'YLim',[0 1.0], 'YTick', 0:0.1:1.0, 'YTickLabel', makeTickLabel(0:0.1:1.0,0.2), 'TickDir', 'out');
xlabel('Motion strength (%coh)');
ylabel('Probability sure target');
h = legend({'nostim','stim'});
set(h,'FontSize',9);
legend('boxoff');
        %Pright
subplot(3,1,2);
hold on;
errorbar(coh_axis/10, pright_strg(:,1), pright_strg_se(:,1), 'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 4);
errorbar(coh_axis/10, pright_nostrg(:,1), pright_nostrg_se(:,1), 'o', 'Color', 'b', 'MarkerFaceColor', 'w', 'MarkerSize', 4);
errorbar(coh_axis/10, pright_strg(:,2), pright_strg_se(:,2), 'o', 'Color', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 4);
errorbar(coh_axis/10, pright_nostrg(:,2), pright_nostrg_se(:,2), 'o', 'Color', 'r', 'MarkerFaceColor', 'w', 'MarkerSize', 4);
plot(g_coh*100, g_pright_strg(:,1), 'b-');
plot(g_coh*100, g_pright_nostrg(:,1), 'b-.');
plot(g_coh*100, g_pright_strg(:,2), 'r-');
plot(g_coh*100, g_pright_nostrg(:,2), 'r-.');
plot([g_coh(1) g_coh(end)]*100,[0.5 0.5], 'k--');
plot([0 0],[0 1],'k--');
set(gca, 'XLim', [-60 60], 'XTick', -60:10:60, 'XTickLabel', makeTickLabel(-60:10:60,20), ...
         'YLim',[-0.02 1.02], 'YTick', 0:0.1:1, 'YTickLabel', makeTickLabel(0:0.1:1,0.2), 'TickDir', 'out');
xlabel('Motion strength (%coh)');
ylabel('Probability right');
h = legend({'w/ stg','w/o strg'},'Location','NorthWest');
set(h,'FontSize',9);
legend('boxoff');

subplot(3,1,3); axis off;
filename2 = filename;
filename2(filename2=='_')= ' ';
if ~exist('modelLL','var')
    modelLL(2) = m2;
end

%         1  2 3     4      5    6  7    8        9      10        11     12     13
      %  [k  B theta stmeff bias dk dsig sig2beta weber guess_rate theta2 dtheta dtheta2] %%%%%%

text(0.1,0.3, sprintf( ... 
    'file: %s\n\t\t\t\t\t\tk= %f\n\t\t\t\t\t\tB= %f\n\t\t\t\t\t\ttheta= %f\n\t\t\t\t\t\tdcohstim= %f\n\t\t\t\t\t\tdcohbias= %f\n\t\t\t\t\t\tdk= %f\n\t\t\t\t\t\tdsigma2= %f\n\t\t\t\t\t\tsigma2beta= %f\n\t\t\t\t\t\tweber= %f\n\t\t\t\t\t\tguess   rate= %f\n\t\t\t\t\t\ttheta2= %f\n\t\t\t\t\t\tdtheta= %f\n\t\t\t\t\t\tdtheta2= %f\n\nerr= %f\nerr2= %f\nnum trials = %d\ndur set = [%d %d]', ... 
    filename2, ...
    modelParam(2).final(1), ...
    modelParam(2).final(2), ...
    modelParam(2).final(3), ...
    modelParam(2).final(4), ...
    modelParam(2).final(5), ...
    modelParam(2).final(6), ...
    modelParam(2).final(7), ...
    modelParam(2).final(8), ...
    modelParam(2).final(9), ...
    modelParam(2).final(10), ...
    modelParam(2).final(11), ...
    modelParam(2).final(12), ...
    modelParam(2).final(13), ...
    -modelLL(2), ...
    err2, ...
    length(dur_(dur_>dur_setQ{q}(1) & dur_<dur_setQ{q}(2))), ...
    dur_setQ{q}(1), ...
    dur_setQ{q}(2)  ));

if savePDF
    print('-dpdf',[folder filename '_diffFitNew_q=' num2str(q)]);
end

% ********** End, make the figure *********

%********** Begin paper figs *********
if paper_figs
% CAREFUL! WILL REPLACE FIGS IN ILLUSTRATOR FILE

% talk figures: larger fonts/symbols, no text reports, etc.
plot_PsureFit = 1;
% PMF
figure; set(gcf, 'Color', [1 1 1], 'Position', [100 100 400 350], 'PaperPositionMode', 'auto'); hold on;

if fourCurves
    hsym1 = plot(g_coh*100, g_pright_nostrg(:,1), 'b--', 'LineWidth', 2);
    hsym2 = plot(g_coh*100, g_pright_nostrg(:,2), 'r--', 'LineWidth', 2);
    hsym3 = plot(g_coh*100, g_pright_strg(:,1), 'b-', 'LineWidth', 2);
    hsym4 = plot(g_coh*100, g_pright_strg(:,2), 'r-', 'LineWidth', 2);    [hsym1,hxe,hye] = errorbar2(coh_axis/10, pright_nostrg(:,1), 0, pright_nostrg_se(:,1), 'o');
    set(hsym1,'Color','b','MarkerSize',11,'MarkerFaceColor','w','MarkerEdgeColor','b');
    set(hye,'LineWidth',1,'Color','b');
    [hsym2,hxe,hye] = errorbar2(coh_axis/10, pright_nostrg(:,2), 0, pright_nostrg_se(:,2), 'o');
    set(hsym2,'MarkerSize',11,'Color','r','MarkerFaceColor','w','MarkerEdgeColor','r');
    set(hye,'LineWidth',1,'Color','r');
    [hsym3,hxe,hye] = errorbar2(coh_axis/10, pright_strg(:,1), 0, pright_strg_se(:,1), 'o');
    set(hsym3,'Color','b','MarkerSize',11,'MarkerFaceColor','b','MarkerEdgeColor','b');
    set(hye,'LineWidth',1,'Color','b');
    [hsym4,hxe,hye] = errorbar2(coh_axis/10, pright_strg(:,2), 0, pright_strg_se(:,2), 'o');
    set(hsym4,'MarkerSize',11,'Color','r','MarkerFaceColor','r','MarkerEdgeColor','r');
    set(hye,'LineWidth',1,'Color','r');
else
    plot(g_coh*100, g_pright_all(:,1), 'b-', 'LineWidth', 2);
    plot(g_coh*100, g_pright_all(:,2), 'r-', 'LineWidth', 2);
    [hsym3,hxe,hye] = errorbar2(coh_axis/10, pright_all(:,1), 0, pright_all_se(:,1), 'o');
    set(hsym3,'Color','b','MarkerSize',11,'MarkerFaceColor','b','MarkerEdgeColor','b');
    set(hye,'LineWidth',1,'Color','b');
    [hsym4,hxe,hye] = errorbar2(coh_axis/10, pright_all(:,2), 0, pright_all_se(:,2), 'o');
    set(hsym4,'MarkerSize',11,'Color','r','MarkerFaceColor','r','MarkerEdgeColor','r');
    set(hye,'LineWidth',1,'Color','r');
end

% % for S3 (or maybe all?):
% plot([0 0],[0 1],'k--','LineWidth',3);
% plot([-53 0],[0.5 0.5],'k--','LineWidth',3);

set(gca, 'XLim', [-53 53], 'XTick', -50:25:50, 'XTickLabel', makeTickLabel(-50:25:50,25), ...
         'YLim',[-0.02 1.02], 'YTick', 0:0.1:1, 'YTickLabel', makeTickLabel(0:0.1:1,0.2), ...
         'TickDir', 'out');
changeAxesFontSize(gca, 22, 22);
if fourCurves
    h=legend([hsym1 hsym3 hsym2 hsym4],'-µS, -T_s','-µS, +T_s','+µS, -T_s','+µS, +T_s','Location', 'Southeast');
%     h=legend([hsym1 hsym3 hsym2 hsym4],'-   coh, -T_s','-   coh, +T_s','+   coh, -T_s','+   coh, +T_s','Location','Southeast');
else
%     h=legend([hsym3 hsym4],'no µS','µS','Location','Northwest');
end
set(h,'FontSize',16);
ylabel('Proportion preferred choices');
xlabel('Motion strength (% coh)');

legend('boxoff'); 

% saveas(gcf,[filename '_diffFit_PMF_wstim_both_q=' num2str(q)],'epsc');
% saveas(gcf,[filename '_diffFit_PMF_wstim_both_q=' num2str(q) '_sigmult'],'epsc');
% saveas(gcf,['FINAL_diffFit_PMF_5param_fit_q=' num2str(q)],'epsc');
% saveas(gcf,['FINAL_diffFit_PMF_bestfit_q=' num2str(q)],'epsc');
% saveas(gcf,['FINAL_diffFit_PMF_bestfit_PSEUDO_q=' num2str(q)],'epsc');
% saveas(gcf,['FINAL_diffFit_PMF_HIGHCURR_q=' num2str(q)],'epsc');
% saveas(gcf,['FINAL_diffFit_PMF_HIGHCURR_nothingisfinal_q=' num2str(q)],'epsc');
% saveas(gcf,['diffFit_PMF_IBIM_q=' num2str(q)],'epsc');
% saveas(gcf,['tempPMF_q=' num2str(q)],'epsc');
% saveas(gcf,['highcurr_best_separated_PMF'],'epsc');

export_fig 'tempPMF' -pdf

% Psure
%         g_pstrg(:,1) = psure_fcn_b6(Psure_beta_indic_full,g_coh,0); % for logistic/gaussian
%         g_pstrg(:,2) = psure_fcn_b6(Psure_beta_indic_full,g_coh,1);
figure; set(gcf, 'Color', [1 1 1], 'Position', [100 100 400 350], 'PaperPositionMode', 'auto'); hold on;
if plot_PsureFit
    h1 = plot(g_coh*100, g_pstrg(:,1), 'b-', 'LineWidth', 2);
    h2 = plot(g_coh*100, g_pstrg(:,2), 'r-', 'LineWidth', 2);
end
[hsym1,hxe,hye] = errorbar2(coh_axis/10, pstrg(:,1), 0, pstrg_se(:,1), 'o');
set(hsym1,'MarkerSize',11,'MarkerFaceColor','b','MarkerEdgeColor','b');
set(hye,'LineWidth',1,'Color','b');
[hsym2,hxe,hye] = errorbar2(coh_axis/10, pstrg(:,2), 0, pstrg_se(:,2), 'o');
set(hsym2,'MarkerSize',11,'MarkerFaceColor','r','MarkerEdgeColor','r');
set(hye,'LineWidth',1,'Color','r');

% for S3:
% plot([0 0],[0 0.8],'k--','LineWidth',2);

% h=legend([hsym1 hsym2 h2],'no µS','µS','prediction','Location','Northeast');
if strcmp(stim_type(1:7),'optstim')
    h=legend([hsym1 hsym2],'no laser','laser on','Location','Northwest');
else
    h=legend([hsym1 hsym2],'no µS','µS','Location','Northeast');
end



% h=legend([hsym1 hsym2],'no   coh','  coh','Location','Northeast');
% h=legend(h2,'Model prediction','Location','Northwest');
set(h,'FontSize',16);

legend('boxoff');

set(gca, 'XLim', [-53 53], 'XTick', -50:25:50, 'XTickLabel', [], ... %makeTickLabel(-50:25:50,25), ...
         'YLim',[0 1], 'YTick', 0:0.1:1, 'YTickLabel', makeTickLabel(0:0.1:1,0.2), ...
         'TickDir', 'out');
% if q==1
    ylabel('Proportion sure-bet choices');
% else
%     set(gca, 'YTickLabel', []);
% end

% temp

% xlabel('Motion strength (% coh)');

changeAxesFontSize(gca, 22, 22);


% saveas(gcf,[filename '_diffFit_Psure_both_wfits_q=' num2str(q)],'epsc');
% saveas(gcf,[filename '_diffFit_Psure_both_wfits_q=' num2str(q) '_sigmult'],'epsc');
% saveas(gcf,['FINAL_diffFit_Psure_5param_pred_q=' num2str(q)],'epsc');
% saveas(gcf,['FINAL_diffFit_Psure_bestfit_q=' num2str(q)],'epsc');
% saveas(gcf,['FINAL_diffFit_Psure_bestfit_PSEUDO_q=' num2str(q)],'epsc');
% saveas(gcf,['FINAL_diffFit_Psure_HIGHCURR_q=' num2str(q)],'epsc');
% saveas(gcf,['diffFit_Psure_HIGHCURR_nothingisfinal_q=' num2str(q)],'epsc');
% saveas(gcf,['diffFit_Psure_IBIM_q=' num2str(q)],'epsc');
% saveas(gcf,['tempPsure_q=' num2str(q)],'epsc');
% saveas(gcf,['highcurr_best_separated_Psure'],'epsc');
% saveas(gcf,['PRED_sig+k_allDurs_Psure'],'epsc');

export_fig 'tempPsure' -pdf


end

end

if saveMAT
    save([folder filename '_diffFit' strategyIndex '.mat']);
end
