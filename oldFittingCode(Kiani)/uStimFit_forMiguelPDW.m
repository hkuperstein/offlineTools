clear all

% if isempty(gcp)
%     cmp = computer;
%     if strcmp(cmp,'MACI64')
%         parpool;
%     elseif strcmp(cmp,'GLNXA64')
%         parpool('metheny')
%     else
%         error('cannot identify computer to set up matlabpool');
%     end
% end

% filename = 'FIRA_uStim_bothmonks_All_N=62'; stim_type = 'elestim';
% load(['[[ENTER FOLDER HERE]]' filename '.mat'], 'FIRA');


clear all
load Hanzo_Confidence

validTrials = ~isnan(signCoherence) & ~isnan(RT) & ~isnan(choice) & ...
              ~isnan(correct) & ~isnan(confidenceChoice);
          
% rename some vars:
D.strength = signCoherence(validTrials);
D.dur = round(RT(validTrials)*1000);
D.choice = choice(validTrials); 
D.correct = correct(validTrials);
D.conf = confidenceChoice(validTrials);

stimtrg = 2;
coh_ = D.strength;
dur_ = D.dur;
    % assume half the trials are strg and generate new choice vector
    % according to conf
strg_ = randsample([0 1], length(coh_), 'true');
cho_trg_ = nan(size(D.choice));
cho_trg_(strg_==1 & D.conf==0) = 3; % low conf = sure bet
cho_trg_(strg_==1 & D.conf==1) = D.choice(strg_==1 & D.conf==1)+1;
cho_trg_(strg_==0) = D.choice(strg_==0)+1;
cor_trg_ = sign(coh_)+1;
cor_trg_(cor_trg_==0) = round(rand)+1;

stim_ = false(size(cho_trg_));
coh_set = unique(coh_);  %how to group coherences
coh_axis = coh_set;      %what coherence is representative of each group
coh_set_unsigned = coh_set(6:end);
coh_axis_unsigned = coh_axis(6:end);

clear D;

fitwhat = {'PrightPS', 'PrightPS'}; 

% Current version has the following 13 parameters:
%        1  2 3     4      5    6  7    8        9      10        11     12     13
%       [k  B theta stmeff bias dk dsig sig2beta weber guess_rate theta2 dtheta dtheta2] %%%%%%
fixed = [0  0 0     1      0    1  1    1        1      1         1      1      1];
% ... but the above '1s' are held fixed


% initial guess
guess = [...
    0.3
    30
    0.6
    0
    0
    0
    0
    0
    0
    3.7201e-44
    0.6
    0
    0]';
        
        D.cor_trg = cor_trg_;
        D.cho_trg = cho_trg_;
        D.coh     = coh_;
        D.dur = dur_;
        D.SureT   = strg_;
        D.ustim = stim_;

            feedback = 1;

            clear options;
            options.fitwhat = fitwhat;
            options.flag = 1;
            options.init_bias_influences_marginals = false;
            options.stim_bias_influences_marginals = false;
            options.use_midline_as_criterion = true;    % set to 1 to use conventional mapping: positive decision variable means right choice and negative means left choice
                                                        % set to 0 to depart from the conventional mapping and instead choose based on the sign of logodds_right 
            options.parfor = true;
            options.plot = false;
            options.coh_set = [];
            options.coh_set_freq = [];
            options.yoke_theta = 1;                          % additional parameters
            options.fitMethod = 'fms';
            
        [modelParam(2),modelLL(2),trialData(2)] = fit_Diffusion_posterior_new4(D, guess, fixed, feedback, options);
        expectedPright{2} = trialData(1).expectedPright;
        expectedPS{2} = trialData(1).expectedPS;
        
        
%         I = abs(trialData(2).coh)==0.512 & strg_==1;
%         pstrg = mean(trialData(2).expectedPS(I))
        

%********** End, fit Pcor and Psure *********


%********** Begin, Make smooth Pright and Psure curves based on model parameters *********

fitParam = modelParam(2).final;
err2 = NaN;

    %generate a set of random stimulus durations.  
dursamples = 1000;      % the actual number of simulated trials will be
                        % dursamples * length(g_coh) * 4
                        % the 4 comes from 2[for_Ts] * 2[for_uStim]

if ~exist('dur_','var')
    dur_ = trialData(2).dur;
end
                        
% can median split to show effect of duration
Q(1) = 0; Q(3) = inf; Q(2) = median(dur_);
dur_setQ{1} = [Q(1) Q(2)];
dur_setQ{2} = [Q(2)+1 Q(3)];
dur_setQ{3} = [0 inf];

% for now just use all trials, dur_setQ{3}
q=3;

rand_dur = randsample(dur_, 2*dursamples, 'true');
g_coh = unique([(-0.6:0.04:0.6)'; cat(2,coh_set)']);

% D is a larger dataset based on resampling original data, 
% used only to generate the smooth curves for display purposes
D.coh       = repmat(g_coh', [length(rand_dur)*4 1]);
D.dur   = repmat(repmat(rand_dur',[1 length(g_coh)]), [4 1]);
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

options.coh_set = unique(trialData(2).coh);
options.coh_set_freq = nan(length(coh_set));        
for c = 1 : length(coh_set)
    s=0;
    options.coh_set_freq(c,1) = sum(trialData(2).coh==options.coh_set(c)&trialData(2).ustim==s)/length(trialData(2).coh);
end

options.plot = 0;
% options.flag = 1;

% call the fitting routine one last time with params fixed
[m1,m2,D] = fit_Diffusion_posterior_new4(D, fitParam, ones(size(fitParam)), 0, options);

g_pright_nostrg = nan(length(g_coh), 1);
g_pright_strg = nan(length(g_coh), 1);
g_pright_all = nan(length(g_coh), 1);
g_pstrg = nan(length(g_coh), 1);
for c = 1 : length(g_coh)
    I = D.coh==g_coh(c) & D.ustim==s & D.SureT==0;
    g_pright_nostrg(c,1) = nanmean(D.expectedPright(I));
    I = D.coh==g_coh(c) & D.ustim==s & D.SureT==1;
    g_pright_strg(c,1) = nanmean(D.expectedPright(I));
    g_pstrg(c,1) = nanmean(D.expectedPS(I));
    I = D.coh==g_coh(c) & D.ustim==s;
    g_pright_all(c,1) = nanmean(D.expectedPright(I));
end

%********** End, Make smooth Pright and Psure curves based on model parameters *********
%%
    
%     %********** re-extract some variables from the data *********
%     all_ind = getTrialIndex_certstim(FIRA,[],[],[],dur_setQ{q},[],[0 1 2]); 
%     stimtrg = 2;
%     coh_set = {-512, -256, -128, -64, [-32 -16], 0, [16 32], 64, 128, 256, 512};  %how to group coherences
%     coh_axis = [-512, -256, -128, -64, -32, 0, 32, 64, 128, 256, 512];      %what coherence is representative of each group
%     coh_ = FIRA{2}(all_ind,IND('dot_coh'));
%     dur_ = FIRA{2}(all_ind,IND('dot_off')) - FIRA{2}(all_ind,IND('dots_on'));       
%     strg_ = ~isnan(FIRA{2}(all_ind,getFIRA_columnByName('s_trg')));
%     cor_ = FIRA{2}(all_ind,IND('correct'));                                         
%     cor_trg_ = FIRA{2}(all_ind,IND('cor_trg'));
%     cho_trg_ = FIRA{2}(all_ind,IND('cho_trg'));
%     coh_(cor_trg_~=stimtrg) = -coh_(cor_trg_~=stimtrg);
%     if ~isempty(FIRA{2}(all_ind,IND('pseudostim'))) && sum(FIRA{2}(all_ind,IND('pseudostim'))) > 5
%         stim_ = FIRA{2}(all_ind,IND('pseudostim'));
%     else
%         stim_ = FIRA{2}(all_ind,IND('elestim'));
%     end
    
        %calculate p(right), p(strg), etc
    pright_strg = nan(length(coh_set), 1);
    pcorr_strg = nan(length(coh_set), 1);
    n_strg = nan(size(pright_strg));
    pright_nostrg = nan(length(coh_set), 1);
    pcorr_nostrg = nan(length(coh_set), 1);
    n_nostrg = nan(size(pright_nostrg));
    pstrg = nan(length(coh_set), 1);
    nstrg = nan(size(pstrg));
    pright_all = nan(length(coh_set), 1);

        %find trials with sure target
    J = strg_==1;

    %calculate Pright and Pstrg for signed coherences 
I = J & ismember(cho_trg_,[1 2]);
[pright_strg(:,1), pright_strg_se(:,1)] = calcGroupMean(cho_trg_(I)==stimtrg, coh_(I), coh_set, 'binary');
I = J;
[pstrg(:,1), pstrg_se(:,1)] = calcGroupMean(cho_trg_(I)==3, coh_(I), coh_set, 'binary');
I = ~J & ismember(cho_trg_,[1 2]);
[pright_nostrg(:,1), pright_nostrg_se(:,1)] = calcGroupMean(cho_trg_(I)==stimtrg, coh_(I), coh_set, 'binary');
I = ismember(cho_trg_,[1 2]);
[pright_all(:,1), pright_all_se(:,1)] = calcGroupMean(cho_trg_(I)==stimtrg, coh_(I), coh_set, 'binary');

%********** plot *********

%    Pright and Psure as a function of signed coherence 
figure(111+q); clf;
set(gcf, 'Color', [1 1 1], 'Position', [300+100*q 300 600 850], 'PaperPositionMode', 'auto');
        %Psure
subplot(3,1,1);
hold on;
errorbar(coh_axis, pstrg(:,1), pstrg_se(:,1), 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
plot(g_coh, g_pstrg(:,1), 'k-');
    plot([g_coh(g_pstrg(:,1)==max(g_pstrg(:,1))) g_coh(g_pstrg(:,1)==max(g_pstrg(:,1)))], [0 max(g_pstrg(:,1))], 'g-', 'LineWidth', 2);
plot([0 0],[0 1],'k--');
set(gca, 'XLim', [-0.6 0.6], 'XTick', -0.6:0.1:0.6, 'XTickLabel', makeTickLabel(-0.6:0.1:0.6,0.2), ...
         'YLim',[0 1.0], 'YTick', 0:0.1:1.0, 'YTickLabel', makeTickLabel(0:0.1:1.0,0.2), 'TickDir', 'out');
xlabel('Motion strength (%coh)');
ylabel('Probability sure target');
h = legend({'nostim','stim'});
set(h,'FontSize',9);
legend('boxoff');
        %Pright
subplot(3,1,2);
hold on;
errorbar(coh_axis, pright_strg(:,1), pright_strg_se(:,1), 'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 4);
errorbar(coh_axis, pright_nostrg(:,1), pright_nostrg_se(:,1), 'o', 'Color', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 4);
plot(g_coh, g_pright_strg(:,1), 'b-');
plot(g_coh, g_pright_nostrg(:,1), 'r-');
plot([g_coh(1) g_coh(end)],[0.5 0.5], 'k--');
plot([0 0],[0 1],'k--');
set(gca, 'XLim', [-0.6 0.6], 'XTick', -0.6:0.1:0.6, 'XTickLabel', makeTickLabel(-0.6:0.1:0.6,0.2), ...
         'YLim',[-0.02 1.02], 'YTick', 0:0.1:1, 'YTickLabel', makeTickLabel(0:0.1:1,0.2), 'TickDir', 'out');
xlabel('Motion strength (%coh)');
ylabel('Probability right');
h = legend({'w/ stg','w/o strg'},'Location','NorthWest');
set(h,'FontSize',9);
legend('boxoff');

subplot(3,1,3); axis off;
if ~exist('modelLL','var')
    modelLL(2) = m2;
end

%         1  2 3     4      5    6  7    8        9      10        11     12     13
      %  [k  B theta stmeff bias dk dsig sig2beta weber guess_rate theta2 dtheta dtheta2] %%%%%%

text(0.1,0.3, sprintf( ... 
    'file: %s\n\t\t\t\t\t\tk= %f\n\t\t\t\t\t\tB= %f\n\t\t\t\t\t\ttheta= %f\n\t\t\t\t\t\tdcoh-stim= %f\n\t\t\t\t\t\tdcoh-bias= %f\n\t\t\t\t\t\tdk= %f\n\t\t\t\t\t\tdsigma= %f\n\t\t\t\t\t\tsigma2-beta= %f\n\t\t\t\t\t\tweber= %f\n\t\t\t\t\t\tguessrate= %f\n\t\t\t\t\t\ttheta2= %f\n\t\t\t\t\t\tdtheta= %f\n\t\t\t\t\t\tdtheta2= %f\n\nerr= %f\nerr2= %f\nnum trials = %d\ndur set = [%d %d]', ... 
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

% ********** End, make the figure *********



