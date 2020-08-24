% function uStimFit(filename, PsMag)

% low-curr, basic model
LLnull = -35601.883515; % 
df = 6;
n = 53134;
BICnull = -2*LLnull + df*log(n)

LLalt = -35540.220345 % 
df = 9;
n = 53134;
BICalt = -2*LLalt + df*log(n)

deltaBIC = BICnull - BICalt
AltSigBetter = (BICalt < BICnull)
pval_LR = chi2pdf(2*(LLalt-LLnull),df)

%%


% LL = -333803.986085; df=7; % kmult and sigmult (same as all three) *** WINNER
% %     load fitParam_4-19-13_allthree_BESTFIT
% % LL = -333814.654082; df=7; % kmult and bmult
% %     load fitParam_4-21-13_bmult_kmult
% % LL = -333821.717633; df=6; % kmult only --> does nothing for Psure height but handles width, and better at capturing PMF slope
% %     load fitParam_4-21-13_kmultONLY
% % LL = -333821.799898; df=7; % sigmult and bmult
% %     load fitParam_4-11-13_bothfreeB;
% % LL = -333821.841382; df=6; % sigmult only --> big eff on Psure height but cannot explain all of slope eff
% %     load fitParam_4-11-13_sigmult_fminsearch
% % LL = -333827.831191; df=7; % sigmult and bmult, bmult forced low
% %     load fitParam_4-21-13_bmult_sigmult_bLow
% % LL = -333881.215040; df=6; % bmult only --> affects mainly Psure height, a tiny bit of slope (?)
%       ??? 
% % LL = -333887.442740; df=5; % 5-param basic model
%       load fitParam_5-17-13_5param_fit
% 
% % BIC / approx. Bayes Factor:
% df = 7;
% n = 53151;
% BIC = -2*LL + df*log(n);
% disp(num2str(BIC));


% % Bmult and sigmaMult free versus both fixed
% LLnull = -333889.465359; LLalt = -333821.799898; df = 2;  
% 
% % Bmult only free versus sigmaMult only free (latter is alt)
% LLnull = -333882.107327; LLalt = -333821.841382; df = 0; %% ONLY BIC, NOT LLR (NOT NESTED)  
% 
% % Kmult only free versus sigmaMult only free (latter is alt)
% LLnull = -333824.5163; LLalt = -333821.841382; df = 0; %% ONLY BIC, NOT LLR (NOT NESTED)  


% % likelihood ratio test:
% pval = chi2pdf(2*(LLalt-LLnull),df)
% 
% 
% % BIC / approx. Bayes Factor:
% dfNull = 6; dfAlt = 6;
% n = 53831;
% 
% BICnull = -2*LLnull + dfNull*log(n);
% BICalt = -2*LLalt + dfAlt*log(n);
% AltSigBetter = (BICalt < BICnull)
% BayesFactor = exp(0.5*(BICnull-BICalt))

clear all
close all

% stim_type = 'pseudostim';
stim_type = 'elestim';

% filename = 'FIRA_uStim_bothmonks_All_N=63'; PsMag = NaN; stim_type = 'elestim';
% filename = 'FIRA_uStim_Ike_All_N=39'; PsMag = NaN;
% filename = 'FIRA_uStim_Damien_All_N=24'; PsMag = NaN;

% filename = 'FIRA_PsStim_All_Sessions_N=42'; PsMag = NaN; stim_type = 'pseudostim'; % g_stimeff = PsMag/100;

% filename = 'FIRA_uStim_Damien_highcurr_N=5'; PsMag = NaN;

filename = 'workspace_NoSlopeReduction_N=48'; PsMag = NaN;

indiv_session = 0;

tic

% Modeling the effect of uStim by assuming uStim and bias act as coherence 
% RK, 2/2011

% clear all
% figure(111); close;
% dbstop if error

savePDF = 1;
saveMAT = 1;
paper_figs = 1; fourCurves = 1;

    %load the data. data.mat contains a FIRA structure that is an aggregate of FIRA
    %structures across different sessions. the FIRA structure of each session is modified
    %so that the preferred direction of stimulated neurons always corresponds to the same target (1 or 2).

load(['/Users/crfetsch/Documents/MATLAB/' filename '.mat'], 'FIRA');
% load(['/home/chris/Documents/MATLAB/' filename '.mat']);
stimtrg=2;

    %extract some of the relevant parameters from FIRA
all_ind = getTrialIndex_certstim(FIRA,[],[],[],[0 inf],[],[0 1 2]);    %all trials
coh_ = FIRA{2}(all_ind,IND('dot_coh'));
dur_ = FIRA{2}(all_ind,IND('dot_off')) - FIRA{2}(all_ind,IND('dots_on'));
delay_ = FIRA{2}(all_ind,IND('fp_off')) - FIRA{2}(all_ind,IND('dot_off'));
strg_delay_ = FIRA{2}(all_ind,IND('s_on')) - FIRA{2}(all_ind,IND('dot_off'));
cor_ = FIRA{2}(all_ind,IND('correct'));
cor_trg_ = FIRA{2}(all_ind,IND('cor_trg'));
cho_trg_ = FIRA{2}(all_ind,IND('cho_trg'));
coh_(cor_trg_~=stimtrg) = -coh_(cor_trg_~=stimtrg);
stim_ = FIRA{2}(all_ind,IND(stim_type));
% session_ = FIRA{1}.ownership;
coh_set = {-512, -256, -128, -64, [-32 -16], 0, [16 32], 64, 128, 256, 512};  %how to group coherences
coh_axis = [-512, -256, -128, -64, -32, 0, 32, 64, 128, 256, 512];      %what coherence is representative of each group
coh_set_unsigned = coh_set(6:end);
coh_axis_unsigned = coh_axis(6:end);
dur_(dur_<60) = 60; % kluge a strange zero-dur trial


if indiv_session
    % try to find best initial guesses by fitting analytic
    % solution for basic diffusion model (a la 2006 chapter)
        %calculate Pright and Pstrg for signed coherences 

    % D = dur_>360; % temp  (>330 or <113)

    J = ~isnan(FIRA{2}(all_ind,getFIRA_columnByName('s_trg')));    %find trials with sure target    
    for s = 0 : 1,
        S = stim_==s;
        I = J & S & ismember(cho_trg_,[1 2]);
        [pright_strg(:,s+1), pright_strg_se(:,s+1)] = calcGroupMean(cho_trg_(I)==stimtrg, coh_(I), coh_set, 'binary');
        I = J & S;
        [pstrg(:,s+1), pstrg_se(:,s+1)] = calcGroupMean(cho_trg_(I)==3, coh_(I), coh_set, 'binary');
        I = ~J & S & ismember(cho_trg_,[1 2]);
        [pright_nostrg(:,s+1), pright_nostrg_se(:,s+1)] = calcGroupMean(cho_trg_(I)==stimtrg, coh_(I), coh_set, 'binary');
    end;

    coh = coh_/1000;
    g_coh = (-0.6:0.02:0.6)';

    PMF_func = @(b,coh) 1./(1+exp(-2*b(1)*b(2)*coh + b(3)));
    % Error to minimize is -LL, where likelihood is the joint probability of
    % obtaining the data set conceived as a list of Bernoulli trials.
    % This is just p1 for rightward choice trials and 1-p1 for leftward.  -after RK (see logistfit_se.m)
    PMF_err_func = @(par,choice1,choice2,coh) -(  sum(log(PMF_func(par,coh(choice1)))) + sum(log(1-PMF_func(par,coh(choice2))))  );

    % I = ~J & stim_==0;
    % [beta1,fval,exitflag,output] = fminsearch(@(x) PMF_err_func(x,cho_trg_(I)==stimtrg,cho_trg_(I)==3-stimtrg,coh(I)), [0.25 15 -0.07]);
    % g_pright_nostrg(:,1) = PMF_func(beta1,g_coh);
    % 
    % I = J & stim_==0;
    % [beta2,fval,exitflag,output] = fminsearch(@(x) PMF_err_func(x,cho_trg_(I)==stimtrg,cho_trg_(I)==3-stimtrg,coh(I)), [0.25 15 -0.07]);
    % g_pright_strg(:,1) = PMF_func(beta2,g_coh);
    % 
    % I = ~J & stim_==1;
    % [beta3,fval,exitflag,output] = fminsearch(@(x) PMF_err_func(x,cho_trg_(I)==stimtrg,cho_trg_(I)==3-stimtrg,coh(I)), [0.25 15 0.14]);
    % g_pright_nostrg(:,2) = PMF_func(beta3,g_coh);
    % 
    % I = J & stim_==1;
    % [beta4,fval,exitflag,output] = fminsearch(@(x) PMF_err_func(x,cho_trg_(I)==stimtrg,cho_trg_(I)==3-stimtrg,coh(I)), [0.25 15 0.14]);
    % g_pright_strg(:,2) = PMF_func(beta4,g_coh);

    % without RT, k and B terms can arbitrarily covary;
    % avoid strange values by bounding with fmincon
    lb = [0.05 5 -inf];
    ub = [0.65 65 inf];

    I = ~J & stim_==0;
    beta1 = fmincon(@(x) PMF_err_func(x,cho_trg_(I)==stimtrg,cho_trg_(I)==3-stimtrg,coh(I)), [0.25 15 1], [], [], [], [], lb, ub);
    g_pright_nostrg(:,1) = PMF_func(beta1,g_coh);

    I = J & stim_==0;
    beta2 = fmincon(@(x) PMF_err_func(x,cho_trg_(I)==stimtrg,cho_trg_(I)==3-stimtrg,coh(I)), [0.25 15 1], [], [], [], [], lb, ub);
    g_pright_strg(:,1) = PMF_func(beta2,g_coh);

    I = ~J & stim_==1;
    beta3 = fmincon(@(x) PMF_err_func(x,cho_trg_(I)==stimtrg,cho_trg_(I)==3-stimtrg,coh(I)), [0.25 15 -1], [], [], [], [], lb, ub);
    g_pright_nostrg(:,2) = PMF_func(beta3,g_coh);

    I = J & stim_==1;
    beta4 = fmincon(@(x) PMF_err_func(x,cho_trg_(I)==stimtrg,cho_trg_(I)==3-stimtrg,coh(I)), [0.25 15 -1], [], [], [], [], lb, ub);
    g_pright_strg(:,2) = PMF_func(beta4,g_coh);

    % beta1
    % beta2
    % beta3
    % beta4

    % guess stimeff and bias from logist_fit parameterization
    I = ~J & ismember(cho_trg_,[1 2]); % sure target not shown
    [fits,~,~,sems] = logist_fit([ones(sum(I),1) stim_(I) (coh(I)).*stim_(I) coh(I) cho_trg_(I)==stimtrg],0,0);
    beta_nostrg_indic(:,1) = fits(3:end); se_nostrg_indic(:,1) = sems(3:end);
    % g_pright_nostrg(:,1) = 1./(1+exp(-(beta_nostrg_indic(1) + beta_nostrg_indic(2)*0 + beta_nostrg_indic(3)*g_coh*0 + beta_nostrg_indic(4)*g_coh))); % no stim
    % g_pright_nostrg(:,2) = 1./(1+exp(-(beta_nostrg_indic(1) + beta_nostrg_indic(2)*1 + beta_nostrg_indic(3)*g_coh*1 + beta_nostrg_indic(4)*g_coh))); % stim

    I = J & ismember(cho_trg_,[1 2]); % sure target shown
    [fits,~,~,sems] = logist_fit([ones(sum(I),1) stim_(I) (coh(I)).*stim_(I) coh(I) cho_trg_(I)==stimtrg],0,0);
    beta_strg_indic(:,1) = fits(3:end); se_strg_indic(:,1) = sems(3:end);
    % g_pright_strg(:,1) = 1./(1+exp(-(beta_strg_indic(1) + beta_strg_indic(2)*0 + beta_strg_indic(3)*g_coh*0 + beta_strg_indic(4)*g_coh))); % no stim
    % g_pright_strg(:,2) = 1./(1+exp(-(beta_strg_indic(1) + beta_strg_indic(2)*1 + beta_strg_indic(3)*g_coh*1 + beta_strg_indic(4)*g_coh))); % stim

    g_k = mean([beta1(1) beta3(1)]);
    g_B = 2*mean([beta1(2) beta3(2)]); % why 2x? I have no idea.
    g_stimeff = mean([beta_strg_indic(2)/beta_strg_indic(4) beta_nostrg_indic(2)/beta_nostrg_indic(4)]);
    if strcmp(stim_type,'pseudostim') && ~isnan(PsMag)
        g_stimeff = PsMag/100;
    end
    g_bias = mean([beta_strg_indic(1)/beta_strg_indic(4) beta_nostrg_indic(1)/beta_nostrg_indic(4)]);
    g_theta = 0.5;

    %********** Begin, explore how well the fits corresponds with the monkey's behavior *********

        %Pright and Psure as a function of signed coherence 
    figure(110); clf;
    set(gcf, 'Color', [1 1 1], 'Position', [100 100 440 300], 'PaperPositionMode', 'auto');
            %Psure
    hold on;
    errorbar(coh_axis/10, pright_strg(:,1), pright_strg_se(:,1), 'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 4);
    errorbar(coh_axis/10, pright_nostrg(:,1), pright_nostrg_se(:,1), 'o', 'Color', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 4);
    errorbar(coh_axis/10, pright_strg(:,2), pright_strg_se(:,2), 'd', 'Color', 'b', 'MarkerFaceColor', 'w', 'MarkerSize', 4);
    errorbar(coh_axis/10, pright_nostrg(:,2), pright_nostrg_se(:,2), 'd', 'Color', 'r', 'MarkerFaceColor', 'w', 'MarkerSize', 4);
    plot(g_coh*100, g_pright_strg(:,1), 'b-');
    plot(g_coh*100, g_pright_nostrg(:,1), 'r-');
    plot(g_coh*100, g_pright_strg(:,2), 'b-.');
    plot(g_coh*100, g_pright_nostrg(:,2), 'r-.');
    set(gca, 'XLim', [-60 60], 'XTick', -60:10:60, 'XTickLabel', makeTickLabel(-60:10:60,20), ...
             'YLim',[-0.02 1.02], 'YTick', 0:0.1:1, 'YTickLabel', makeTickLabel(0:0.1:1,0.2), 'TickDir', 'out');
    xlabel('Motion strength (%coh)');
    ylabel('Probability right');
    h = legend({'w/ stg','w/o strg'},'Location','NorthWest');
    set(h,'FontSize',9);
    legend('boxoff');
end

%********** Begin, fit Pcor and Psure *********

% The model has the following free parameters
%   K               the coefficient that converts the motion strength to evidence.
%   B               bound height.
%   logodds_crit    the criterion on log-odds correct that governs choosing the sure target. 
%   weber_fraction  potentially useful, but turns out to be ineffective here. 
%   uStim_eff       the size of the micto-stim effect expressed in effective coherence.
%   init_bias       the size of bias in non-stimulated trials. for fitting purposes it
%                   may be better to get the value from the data and keep this parameter
%                   fixed.
%   thetaOffset     offset applied to logodds_crit ("theta") for ustim trials, assuming (if nonzero)
%                   that ustim feels 'peculilar' (i.e., to explain mismatched Psure and PMF shifts)


    %there are several possible strategies to approach the fitting. define which strategy
    %should be used.
    
%fit Pright and Psure of non-stimulated trials and Pright of stimulated trials to calculate k, B,
%logodds_crit and init_bias, establishing a *prediction* for Psure of stim trials
FitStrategy = 'OneStage-PrightPS,Pright'; strategyIndex = 'A';

% fit Pright and Psure of all trials to calculate k, B, logodds_crit and init_bias
% FitStrategy = 'OneStage-PrightPS,PrightPS'; strategyIndex = 'B';


    %run the fitting procedure for the specified strategy
switch FitStrategy,
    case {'TwoStage-PrightPS,Pright','TwoStage-PrightPS,PrightPS'},
        if strcmp(FitStrategy,'TwoStage-PrightPS,Pright'),
            fitwhat = {'PrightPS', 'Pright'}; 
        else
            fitwhat = {'PrightPS', 'PrightPS'}; 
        end;
            %first stage: fit Pright and Psure of non-stimulated trials to calculate k, B, logodds_crit and init_bias
        guess = [g_k g_B g_theta 0 0 g_bias 0];
        fixed = [0   0   0       1 1 0      1];
        L = find(stim_==0);
        D.CorrTrg    = cor_trg_(L);
        D.Resp    = cho_trg_(L);
        D.Coh     = coh_(L)/1000;      %must be signed coherence, this is different from the Science paper 
        D.StimDur = dur_(L);
        D.Delay   = delay_(L);
        D.SureT   = ~isnan(strg_delay_(L));
        D.uStim = stim_(L);
        [modelParam(1),modelLL(1),trialData(1)] = fit_Diffusion_posterior_PriorAsCoh(D, guess, fixed, fitwhat, [], [], filename);
        expectedPright{1} = trialData(1).expectedPright;
        expectedPS{1} = trialData(1).expectedPS;

            %second stage: fit stimulated trials to calculate uStim_eff and thetaOffset
        guess = modelParam.final; % guess for k, B, theta, bias
        guess(5) = guess_stimeff;  % guess for uStim_eff
        fixed = [1 1 1 1 0 1 1];
        L = find(stim_==1);
        D.CorrTrg    = cor_trg_(L);
        D.Resp    = cho_trg_(L);
        D.Coh     = coh_(L)/1000;      %must be signed coherence, this is different from the Science paper 
        D.StimDur = dur_(L);
        D.Delay   = delay_(L);
        D.SureT   = ~isnan(strg_delay_(L));
        D.uStim = stim_(L);
        [modelParam(2),modelLL(2),trialData(2)] = fit_Diffusion_posterior_PriorAsCoh(D, guess, fixed, fitwhat, [], [], filename);
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
        %%%%%%  [k   B  theta w stmeff bias TO scb kmult Bmult varOffset] %%%%%%
        
%                 damien FITS: 8P(k+B+VO) = 13440.77748 (but VO basically zero; this is via k and b)
%                              7P(km+VO)  = 13434.5459 - BEST! (given df)
%                              7P(km+TO)  = 13444.9374
%                              9P(w/TO)   = 13433.8714

%                    ike FITS: 7P(km+VO) = 320127.6797 - BEST
%                              7P(km+TO) = 320129.513096


%         load fitParam_6-25-13_HIGHCURR_k+sigmult+oldTO; %1599.8591
%         load fitParam_6-25-13_HIGHCURR_k+sigmult+oldTO#2 % 1600.04
        % WHY DOESN'T IBIM DO THE SAME THING  AS OLD-TO? 
        
%         % drastic change incoming:
%         % replace B with sigma as free param, add param(12) which is
%         % coh-dependent variance term, and then stim can modify both the
%         % coh-var term (13) and an additional additive variance offset (11)
%         guess(2) = 1; % sigma -- NEXT: CHANGE B TO 1 AND SIGMA TO ???
%         guess(12) = 1; % coh-var term
%         guess(13) = 0; % coh-var offset
%         guess(11) = 0.25; % var offset

%         load fitParam_6-27-13_SigFree_CohVar_FIT;
%         load fitParam_6-30-13_SigFree_CohVar+km_FIT; guess=fitParam; % 333772.9050
        load fitParam_6-30-13_SigFree_CohVar+km_PRED; guess=fitParam;
%                     PRED, err = 175373.2979 **** BEST PRED

%         load fitParam_6-24-13_7P(k+TO)_FIT; guess=fitParam; % give TO a chance here...
%         guess(2) = 0.9415; guess(12) = 1.0063; guess(13) = 0;
            % and of course it does better (333771.9672) -- CRAP!

        % now we have to see which model gives best fit to the fitted portion under pred method
%         load fitParam_6-30-13_SigFree_TO+km_FIT; guess=fitParam;
            % #1: 175381.7323 (unc and guess(7)=0)
            % #2: 175379.6509 (search and guess(7)=0)
            % #3: 175379.6674 (search and guess(7)~=0)
            % #4: 175381.8975 (unc and guess(7)~=0)
%                 NONE ARE BETTER THAN VO MODEL
%                     whew, at least fitted portion of pred isn't better w. TO
%                     BUT DOES VO ACTUALLY BUDGE IF YOU SET DK LOWER TO START?
%                         (pred w/ favored model)
%                             --> no, but that doesn't matter; fit is worse (175377.9901)

%         load fitParam_7-05-13_SigFree_VO+km_PSEUDO_PRED; guess = fitParam;
        % figure out whether CohVar really helps w pseudo
            % fitParam_7-05-13_SigFree_VO+km_PSEUDO_FIT_CVguess75 = 552504.3496
            % fitParam_7-05-13_SigFree_VO+km_PSEUDO_FIT_CVguess0  = 552500.5539 % answer is NO
            
            
        % to do:
            % -- see if site selection in damien changes the story
%                 DONE: Psure height effect still clearly there for posOnly (as is mismatch)
            % -- fit high-curr (pred) w/ only dSigma to see how bad Psure misses
                % ... first fitting newest model to highcurr to get a baseline
%                 load fitParam_7-08-13_HIGHCURR_SigFree_9P;% not as good as oldTO, but much better err2
%                 DONE, saved
            % -- CONFIRM fig. s4 still true w/ new bias outside parens

            % -- re-fit Pseudo w/ new model
                % DONE: seems the same, misses a bit, and doesn't do better w. coh-var...
            
            % -- fit each w/ v0 bias
%                   LOW: fitParam_7-09-13_km+VO+x0_noCohVar_FIT: 333776.4373; better! though not by eye
%                   HIGH: 
%                 load fitParam_5-25-13_HIGHCURR_k+sigmult; guess=fitParam;
%                 guess(6)=guess(6)*guess(1); guess(11)=0.03; guess(12)=0; guess(13)=0; guess(14)=0;
                % didn't help asymmetry
%                 load fitParam_temp;
%                 fitParam(14) = -0.25;
%                 fitParam(6) = 0.04; guess=fitParam;
                    % this starting point tries to get closer to assym,
                    % but fit moves each bias term in the opposite dir
            
            % -- fit NoSlopeReduction (7P model)
%                 load fitParam_7-11-13_NoSlopeReduction_7P
                % first dk was large! err = 227912.7526
                % try new starting pt:
%                 fitParam(11) = 0.11; fitParam(9) = 0.9;
%                 guess=fitParam; % BETTER!
                    
%         1  2    3   4    5    6   7   8  9  10 11   12     13 14
      %  [k sig theta w stmeff bias TO scb km Bm VO coh-var CVO x0] %%%%%%        
 fixed = [0  0    0   1    0    0   1   1  0  1  0    1      1  1]; % 1  1]; % last two are for bound change (can be omitted)
        
        if indiv_session
            guess = [g_k g_B g_theta 0 g_stimeff g_bias 0 0 1 1 1];
        end
        D.CorrTrg    = cor_trg_;
        D.Resp    = cho_trg_;
        D.Coh     = coh_/1000;      %must be signed coherence, this is different from the Science paper 
        D.StimDur = dur_;
        D.Delay   = delay_;
        D.SureT   = ~isnan(strg_delay_);
        D.uStim = stim_;
        [modelParam(2),modelLL(2),trialData(2)] = fit_Diffusion_posterior_PriorAsCoh(D, guess, fixed, fitwhat, [], [], filename);
        expectedPright{2} = trialData(1).expectedPright;
        expectedPS{2} = trialData(1).expectedPS;
end;

%********** End, fit Pcor and Psure *********


%********** Begin, Make smooth Pright and Psure curves based on model parameters *********

fitParam = modelParam(2).final;
global err2
    
%% temp: hand-tune params -- comment out this line (cell begin) to "HERE RESUMES NORMAL CODE", below, to resume normal operation

clear all
close all

filename = 'FIRA_uStim_bothmonks_All_N=63'; PsMag = NaN;
% filename = 'FIRA_PsStim_All_Sessions_N=42'; PsMag = NaN;
% filename = 'FIRA_uStim_Damien_highcurr_N=5'; PsMag = NaN;
% filename = 'FIRA_NoSlopeReduction_N=48'; PsMag = NaN;

paper_figs = 1;
fourCurves = 1;

% load fitParam_7-11-13_NoSlopeReduction_7P

load fitParam_6-30-13_SigFree_CohVar+km_FIT
fitParam(2)=30; fitParam(14)=0;


% see codebank(41) for junk that was here


modelParam(2).final = fitParam;
if numel(fitParam)==11
    modelParam(2).final(12)=NaN;
    modelParam(2).final(13)=NaN;
end

savePDF = 1;
saveMAT = 0;

load(['/Users/crfetsch/Documents/MATLAB/' filename '.mat'], 'FIRA');
stimtrg=2;

if isnan(PsMag)
    stim_type = 'elestim'; guess_stimeff = 0.2;
else
    stim_type = 'pseudostim'; guess_stimeff = PsMag/100;
end

    %extract some of the relevant parameters from FIRA
all_ind = getTrialIndex_certstim(FIRA,[],[],[],[0 inf],[],[0 1 2]);    %all trials
coh_ = FIRA{2}(all_ind,IND('dot_coh'));
dur_ = FIRA{2}(all_ind,IND('dot_off')) - FIRA{2}(all_ind,IND('dots_on'));
delay_ = FIRA{2}(all_ind,IND('fp_off')) - FIRA{2}(all_ind,IND('dot_off'));
strg_delay_ = FIRA{2}(all_ind,IND('s_on')) - FIRA{2}(all_ind,IND('dot_off'));
cor_ = FIRA{2}(all_ind,IND('correct'));
cor_trg_ = FIRA{2}(all_ind,IND('cor_trg'));
cho_trg_ = FIRA{2}(all_ind,IND('cho_trg'));
coh_(cor_trg_~=stimtrg) = -coh_(cor_trg_~=stimtrg);
stim_ = FIRA{2}(all_ind,IND(stim_type));
% session_ = FIRA{1}.ownership;
coh_set = {-512, -256, -128, -64, [-32 -16], 0, [16 32], 64, 128, 256, 512};  %how to group coherences
coh_axis = [-512, -256, -128, -64, -32, 0, 32, 64, 128, 256, 512];      %what coherence is representative of each group
coh_set_unsigned = coh_set(6:end);
coh_axis_unsigned = coh_axis(6:end);

dur_(dur_<60) = 60; % kluge a strange zero-dur trial

trialData(2).StimDur = dur_;
fitwhat = {'PrightPS', 'PrightPS'}; strategyIndex = 'B';
% fitwhat = {'PrightPS', 'Pright'}; strategyIndex = 'A';

% first recompute the likelihood of the real data given these params
D.CorrTrg    = cor_trg_;
D.Resp    = cho_trg_;
D.Coh     = coh_/1000;      %must be signed coherence, this is different from the Science paper 
D.StimDur = dur_;
D.Delay   = delay_;
D.SureT   = ~isnan(strg_delay_);
D.uStim = stim_;
[~,LL,~] = fit_Diffusion_posterior_PriorAsCoh(D, fitParam, ones(size(fitParam)), fitwhat, [], [], filename);
modelLL(2) = LL;
global err2
disp(['err2 = ' num2str(err2)]);
pause(0);

% was debugging RK vs. me here; see codebank(38)

%%%%%
%%% HERE RESUMES NORMAL CODE -- comment out above, from here to cell-begin, for normal operation
%%%%%


    %generate a set of random stimulus durations.  
dursamples = 500;       % the actual number of simulated trials will be
                        % dursamples * length(g_coh) * 4
                        % the 4 comes from 2[for_Ts] * 2[for_uStim]

% Q(1) = 0; Q(5) = inf; Q(2:4) = quantile(trialData(2).StimDur,0.25:0.2:0.75); % quartiles
% dur_setQ{1} = [Q(1) Q(2)];
% dur_setQ{2} = [Q(2) Q(3)];
% dur_setQ{3} = [Q(3) Q(4)];
% dur_setQ{4} = [Q(4) Q(5)];
% dur_setQ{5} = [0 inf];

% OR: median split
Q(1) = 0; Q(3) = inf; Q(2) = median(trialData(2).StimDur);
dur_setQ{1} = [Q(1) Q(2)];
dur_setQ{2} = [Q(2)+1 Q(3)];
dur_setQ{3} = [0 inf];

% for q = 1:5 % run/plot for 4 dur quantiles as well as full dataset
for q = 3 % run/plot for 2 halves as well as full dataset

if q==3 %5
    rand_dur = randsample(trialData(2).StimDur, 2*dursamples, 'true');
else
    rand_dur = randsample(trialData(2).StimDur(trialData(2).StimDur>Q(q) & trialData(2).StimDur<Q(q+1)), 2*dursamples, 'true');
end

g_coh = unique([(-0.6:0.04:0.6)'; cat(2,coh_set{:})'/1000]);

D.Coh       = repmat(g_coh', [length(rand_dur)*4 1]);
D.StimDur   = repmat(repmat(rand_dur',[1 length(g_coh)]), [4 1]);
D.CorrTrg      = ones(size(D.Coh));
D.CorrTrg(D.Coh<0) = 2;
D.Resp      = D.CorrTrg;
D.Delay     = repmat(1200, size(D.Coh));
D.SureT     = ones(size(D.Coh));
D.SureT(1:2*dursamples,:) = 0;
D.uStim = ones(size(D.Coh));
D.uStim([1:dursamples, 2*dursamples+(1:dursamples)],:) = 0;
    %make them vertical vectors
D.Coh = D.Coh(:);
D.StimDur = D.StimDur(:);
D.CorrTrg = D.CorrTrg(:);
D.Resp = D.Resp(:);
D.Delay = D.Delay(:);
D.SureT = D.SureT(:);
D.uStim = D.uStim(:);

orig_coh_set = [-0.512 -0.256 -0.128 -0.064 -0.032 0 0.032 0.064 0.128 0.256 0.512];
orig_coh_freq = [1/12*ones(1,5), 1/12*2, 1/12*ones(1,5)];

[m1,m2,D] = fit_Diffusion_posterior_PriorAsCoh(D, fitParam, ones(size(fitParam)), fitwhat, orig_coh_set, orig_coh_freq, filename);

g_pright_nostrg = nan(length(g_coh), 2);
g_pright_strg = nan(length(g_coh), 2);
g_pright_all = nan(length(g_coh), 2);
g_pstrg = nan(length(g_coh), 2);
for s = 0 : 1
    for c = 1 : length(g_coh)
        I = D.Coh==g_coh(c) & D.uStim==s & D.SureT==0;
        g_pright_nostrg(c,s+1) = nanmean(D.expectedPright(I));
        I = D.Coh==g_coh(c) & D.uStim==s & D.SureT==1;
        g_pright_strg(c,s+1) = nanmean(D.expectedPright(I));
        g_pstrg(c,s+1) = nanmean(D.expectedPS(I));
        I = D.Coh==g_coh(c) & D.uStim==s;
        g_pright_all(c,s+1) = nanmean(D.expectedPright(I));
    end
end

% save diffusion_fit_posterior.mat ...
%         FitStrategy fitwhat modelParam modelLL trialData expectedPright expectedPS ...
%         coh_set coh_axis dur_ rand_dur g_coh g_pright_nostrg g_pright_strg g_pstrg;

% save(['/Users/crfetsch/Documents/MATLAB/' filename '_diffFitOrig.mat'], 'FitStrategy', ...
% 'fitwhat', 'modelParam', 'modelLL', 'trialData', 'expectedPright', 'expectedPS', 'coh_set', ...
% 'coh_axis', 'dur_', 'rand_dur', 'g_coh', 'g_pright_nostrg', 'g_pright_strg', 'g_pstrg');
% 
% if saveMAT
%     save(['/Users/crfetsch/Documents/MATLAB/' filename '_diffFit' strategyIndex num2str(length(fitParam)) '_q=' num2str(q) '.mat']);
% end
%********** End, Make smooth Pright and Psure curves based on model parameters *********

    
    %********** Begin, measure Pcor and Psure *********
    all_ind = getTrialIndex_certstim(FIRA,[],[],[],dur_setQ{q},[],[0 1 2]); 

    coh_ = FIRA{2}(all_ind,IND('dot_coh'));
    dur_ = FIRA{2}(all_ind,IND('dot_off')) - FIRA{2}(all_ind,IND('dots_on'));       
    strg_delay_ = FIRA{2}(all_ind,IND('s_on')) - FIRA{2}(all_ind,IND('dot_off'));
    cor_ = FIRA{2}(all_ind,IND('correct'));                                         
    cor_trg_ = FIRA{2}(all_ind,IND('cor_trg'));
    cho_trg_ = FIRA{2}(all_ind,IND('cho_trg'));
    coh_(cor_trg_~=stimtrg) = -coh_(cor_trg_~=stimtrg);
    if ~isempty(FIRA{2}(all_ind,IND('pseudostim'))) && sum(FIRA{2}(all_ind,IND('pseudostim'))) > 5
        stim_ = FIRA{2}(all_ind,IND('pseudostim'));
    else
        stim_ = FIRA{2}(all_ind,IND('elestim'));
    end
    
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
errorbar(coh_axis/10, pright_nostrg(:,1), pright_nostrg_se(:,1), 'o', 'Color', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 4);
errorbar(coh_axis/10, pright_strg(:,2), pright_strg_se(:,2), 'd', 'Color', 'b', 'MarkerFaceColor', 'w', 'MarkerSize', 4);
errorbar(coh_axis/10, pright_nostrg(:,2), pright_nostrg_se(:,2), 'd', 'Color', 'r', 'MarkerFaceColor', 'w', 'MarkerSize', 4);
plot(g_coh*100, g_pright_strg(:,1), 'b-');
plot(g_coh*100, g_pright_nostrg(:,1), 'r-');
plot(g_coh*100, g_pright_strg(:,2), 'b-.');
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
text(0.1,0.3, sprintf( ... 
    'file: %s\n\t\t\t\t\t\tk= %f\n\t\t\t\t\t\tB= %f\n\t\t\t\t\t\tlogodds crit= %f\n\t\t\t\t\t\tweber fraction= %f\n\t\t\t\t\t\tuStim effect= %f\n\t\t\t\t\t\tinit bias= %f\n\t\t\t\t\t\tthetaOffset= %f\n\t\t\t\t\t\tstimChoiceBias= %f\n\t\t\t\t\t\tstimKmult= %f\n\t\t\t\t\t\tstimBmult= %f\n\t\t\t\t\t\tstimVarOffset= %f\n\t\t\t\t\t\tcoh-var= %f\n\t\t\t\t\t\tCVO= %f\n\t\t\t\t\t\tx0= %f\n\nerr= %f\nerr2= %f\nnum trials = %d\ndur set = [%d %d]', ... 
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
    modelParam(2).final(14), ...
    -modelLL(2), ...
    err2, ...
    length(trialData(2).StimDur(trialData(2).StimDur>dur_setQ{q}(1) & trialData(2).StimDur<dur_setQ{q}(2))), ...
    dur_setQ{q}(1), ...
    dur_setQ{q}(2)  ));

if savePDF
    print('-dpdf',[filename '_diffFitNew' strategyIndex '_q=' num2str(q)]);
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
    h=legend([hsym3 hsym4],'no µS','µS','Location','Northwest');
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
saveas(gcf,['tempPMF_q=' num2str(q)],'epsc');
% saveas(gcf,['highcurr_best_separated_PMF'],'epsc');


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
h=legend([hsym1 hsym2],'no µS','µS','Location','Northeast');


% h=legend([hsym1 hsym2],'no   coh','  coh','Location','Northeast');
% h=legend(h2,'Model prediction','Location','Northwest');
set(h,'FontSize',16); legend('boxoff');

set(gca, 'XLim', [-53 53], 'XTick', -50:25:50, 'XTickLabel', [], ... %makeTickLabel(-50:25:50,25), ...
         'YLim',[0 0.8], 'YTick', 0:0.1:0.8, 'YTickLabel', makeTickLabel(0:0.1:0.8,0.2), ...
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
saveas(gcf,['tempPsure_q=' num2str(q)],'epsc');
% saveas(gcf,['highcurr_best_separated_Psure'],'epsc');
% saveas(gcf,['PRED_sig+k_allDurs_Psure'],'epsc');

end

end

if saveMAT
    save(['/Users/crfetsch/Documents/MATLAB/' filename '_diffFit' strategyIndex '.mat']);
end






% % % 
% % % 
% % % 
% % % 
% % % plot_errorbars = 0;
% % % % build option 1: nostrg nostim -> strg nostim -> both stims
% % % plot(g_coh*100, g_pright_nostrg(:,1), 'r-', 'LineWidth', 3);
% % % if plot_errorbars
% % %     [hsym,hxe,hye] = errorbar2(coh_axis/10, pright_nostrg(:,1), 0, pright_nostrg_se(:,1));
% % %     set(hsym,'Color','b','MarkerSize',14,'MarkerFaceColor','w');
% % %     set(hye,'LineWidth',1,'Color','r');
% % % else
% % %     plot(coh_axis/10, pright_nostrg(:,1), 'or', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
% % % %     title(['trials per point (red): ' num2str(min(n_nostrg(:,1))) '-' num2str(max(n_nostrg(:,1)))],'FontSize',12);
% % % end
% % % set(gca, 'XLim', [-53 53], 'XTick', -50:25:50, 'XTickLabel', makeTickLabel(-50:25:50,25), ...
% % %          'YLim',[-0.02 1.02], 'YTick', 0:0.1:1, 'YTickLabel', makeTickLabel(0:0.1:1,0.2), 'TickDir', 'out');
% % % changeAxesFontSize(gca, 18, 18);
% % % % saveas(gcf,[file '_PMF_nostim_nostrg'],'fig');
% % % % saveas(gcf,[file '_PMF_nostim_nostrg'],'epsc');
% % % 
% % % plot(g_coh*100, g_pright_strg(:,1), 'b-', 'LineWidth', 3);
% % % if plot_errorbars
% % %     [hsym,hxe,hye] = errorbar2(coh_axis/10, pright_strg(:,1), 0, pright_strg_se(:,1));
% % %     set(hsym,'Color','b','MarkerSize',14,'MarkerFaceColor','b');
% % %     set(hye,'LineWidth',1,'Color','b');
% % % else
% % %     plot(coh_axis/10, pright_strg(:,1), 'ob', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
% % % %     title(['trials per point (blue): ' num2str(min(n_strg(:,1))) '-' num2str(max(n_strg(:,1)))],'FontSize',12);
% % % end
% % % % saveas(gcf,[file '_PMF_nostim_both'],'fig');
% % % % saveas(gcf,[file '_PMF_nostim_both'],'epsc');
% % % 
% % % plot(g_coh*100, g_pright_nostrg(:,2), 'r--', 'LineWidth', 3);
% % % plot(g_coh*100, g_pright_strg(:,2), 'b--', 'LineWidth', 3);
% % % if plot_errorbars
% % %     [hsym,hxe,hye] = errorbar2(coh_axis/10, pright_nostrg(:,2), 0, pright_nostrg_se(:,2));
% % %     set(hsym,'MarkerSize',14,'Color','r','MarkerFaceColor','w');
% % %     set(hye,'LineWidth',1,'Color','r');
% % %     [hsym,hxe,hye] = errorbar2(coh_axis/10, pright_strg(:,2), 0, pright_strg_se(:,2));
% % %     set(hsym,'MarkerSize',14,'Color','b','MarkerFaceColor','w');
% % %     set(hye,'LineWidth',1,'Color','b');
% % % else
% % %     plot(coh_axis/10, pright_nostrg(:,2), 'or', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
% % %     plot(coh_axis/10, pright_strg(:,2), 'ob', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
% % % %     title(['trials per point (dashed): ' num2str(min([n_strg(:,2);n_nostrg(:,2)])) '-' num2str(max([n_strg(:,2);n_nostrg(:,2)]))],'FontSize',12);
% % % end
% % % % saveas(gcf,[file '_PMF_wstim_both'],'fig');
% % % 

% 
% 
% sprint_txt = ['%s\t'];
% for o = 1:250 % length just needs to be long enough to cover fields
%     sprint_txt = [sprint_txt '%4.3f\t'];
% end
% outfile = ['/Users/crfetsch/Documents/MATLAB/uStimFit_loop.txt'];
% createfile = 0;
% if (exist(outfile, 'file') == 0)    %file does not yet exist
%     createfile = 1;
% end
% fid = fopen(outfile, 'a');
% if (createfile)
%     fprintf(fid, 'filename\t K\t B\t theta\t weber\t uStim_eff\t init_bias\t thetaOffset\t scbias\t skmult\t modelLL\t PsMag\t');
%     fprintf(fid, '\r\n');
% end
% buff = sprintf(sprint_txt, filename, modelParam(2).final, modelLL(2), PsMag);
% fprintf(fid, '%s', buff);
% fprintf(fid, '\r\n');
% fclose(fid);
% 
toc
