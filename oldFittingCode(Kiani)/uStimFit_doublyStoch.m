clear all
close all

cmp = computer;
if strcmp(cmp,'MACI64')
    folder = '/Users/crfetsch/Documents/MATLAB/';
    if matlabpool('size') == 0
        matlabpool open local4
    end
elseif strcmp(cmp,'GLNXA64')
    folder = '/home/chris/MATLAB/';
    if matlabpool('size') == 0
        parpool('metheny',6);
    end
else
    error('cannot identify computer');
end


filename = 'FIRA_Damien_doublyStoch_N=6';

savePDF = 1;
saveMAT = 1;

load([folder filename '.mat'], 'FIRA');

stimtrg=2;

    %extract some of the relevant parameters from FIRA
all_ind = getTrialIndex_certstim(FIRA,[],[],[],[0 inf],[],[0 1 2]);    %all trials
coh_ = FIRA{2}(all_ind,IND('dot_coh'));
coh_std = FIRA{2}(all_ind,IND('coh_std'));
dur_ = FIRA{2}(all_ind,IND('dot_off')) - FIRA{2}(all_ind,IND('dots_on'));
delay_ = FIRA{2}(all_ind,IND('fp_off')) - FIRA{2}(all_ind,IND('dot_off'));
strg_ = ~isnan(FIRA{2}(all_ind,getFIRA_columnByName('s_trg')));
cor_ = FIRA{2}(all_ind,IND('correct'));
cor_trg_ = FIRA{2}(all_ind,IND('cor_trg'));
cho_trg_ = FIRA{2}(all_ind,IND('cho_trg'));
coh_(cor_trg_~=stimtrg) = -coh_(cor_trg_~=stimtrg);
% % % stim_ = FIRA{2}(all_ind,IND(stim_type));
% session_ = FIRA{1}.ownership;
coh_set = {-512, -256, -128, -64, [-32 -16], 0, [16 32], 64, 128, 256, 512};  %how to group coherences
coh_axis = [-512, -256, -128, -64, -32, 0, 32, 64, 128, 256, 512];      %what coherence is representative of each group
coh_set_unsigned = coh_set(6:end);
coh_axis_unsigned = coh_axis(6:end);
dur_(dur_<60) = 60; % kluge a strange zero-dur trial


% shortcut for now is to treat high-coh_std as 'stim' and try to explain results with just dsigma2
stim_ = coh_std==256;
unique(coh_std)


% what to fit: {on stim trials, on nonstim trials}; the rest will be a prediction
% fitwhat = {'PrightPS', 'PrightPS'}; 
fitwhat = {'PrightPS', 'Pright'}; 

%        1  2 3     4      5    6  7    8        9      10        11     12     13
%       [k  B theta stmeff bias dk dsig2 sig2beta weber guess_rate theta2 dtheta dtheta2] %%%%%%
fixed = [0  0 0     1      0    1  0     0        1     1          1      1      1];

k= 0.47
B= 22
theta= 0.8
dcohstim= 0
dcohbias= 0
dk= 0
dsigma2= 0.8
sigma2beta= 1.3
weber= 0
guessrate= 0
theta2= theta
dtheta= 0
dtheta2= 0

guess = [k B theta dcohstim dcohbias dk dsigma2 sigma2beta weber guessrate theta2 dtheta dtheta2];

%     fixed(:)=1;
        
allfixed = all(fixed==1);

D.cor_trg = cor_trg_;
D.cho_trg = cho_trg_;
D.coh     = coh_/1000;
D.dur = dur_;
D.SureT   = strg_;
D.ustim = stim_;

    clear options;
    options.fitwhat = fitwhat;
    options.yoke_theta = 1; % yoke thetas and dthetas (formerly "flag")
    options.fitMethod = 'fms';
[modelParam(2),modelLL(2),trialData(2)] = fit_Diffusion_posterior_new4(D, guess, fixed, 0, options);
expectedPright{2} = trialData(1).expectedPright;
expectedPS{2} = trialData(1).expectedPS;

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
dursamples = 1000;      % the actual number of simulated trials will be
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
for q = 3 % 1:3 % run/plot for 2 halves as well as full dataset

if q==3 %5
    rand_dur = randsample(dur_, 2*dursamples, 'true');
else
    rand_dur = randsample(dur_(dur_>Q(q) & dur_<Q(q+1)), 2*dursamples, 'true');
end

g_coh = unique([(-0.6:0.04:0.6)'; cat(2,coh_set{:})'/1000]);

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

% options.coh_set = [-0.512 -0.256 -0.128 -0.064 -0.032 0 0.032 0.064 0.128 0.256 0.512]';
% options.coh_set_freq = [1/12*ones(1,5), 1/12*2, 1/12*ones(1,5)]';
% options.coh_set_freq = [options.coh_set_freq options.coh_set_freq];

options.coh_set = unique(trialData(2).coh);
options.coh_set_freq = nan(length(coh_set), 2);        
for s = 0 : 1
    for c = 1 : length(coh_set)
        options.coh_set_freq(c,s+1) = sum(trialData(2).coh==options.coh_set(c)&trialData(2).ustim==s)/length(trialData(2).coh);
    end
end

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
    coh_std = FIRA{2}(all_ind,IND('coh_std'));
    dur_ = FIRA{2}(all_ind,IND('dot_off')) - FIRA{2}(all_ind,IND('dots_on'));       
%     strg_delay_ = FIRA{2}(all_ind,IND('s_on')) - FIRA{2}(all_ind,IND('dot_off'));
        strg_ = ~isnan(FIRA{2}(all_ind,getFIRA_columnByName('s_trg')));
    cor_ = FIRA{2}(all_ind,IND('correct'));                                         
    cor_trg_ = FIRA{2}(all_ind,IND('cor_trg'));
    cho_trg_ = FIRA{2}(all_ind,IND('cho_trg'));
    coh_(cor_trg_~=stimtrg) = -coh_(cor_trg_~=stimtrg);
%     if ~isempty(FIRA{2}(all_ind,IND('pseudostim'))) && sum(FIRA{2}(all_ind,IND('pseudostim'))) > 5
%         stim_ = FIRA{2}(all_ind,IND('pseudostim'));
%     else
%         stim_ = FIRA{2}(all_ind,IND('elestim'));
%     end
    % shortcut for now is to treat high-coh_std as 'stim' and try to explain results with just dsigma2
    stim_ = coh_std==256;
    unique(coh_std)

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
subplot(3,2,[1 2]);
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
subplot(3,2,3); hold on;
errorbar(coh_axis/10, pright_strg(:,1), pright_strg_se(:,1), 'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 4);
errorbar(coh_axis/10, pright_nostrg(:,1), pright_nostrg_se(:,1), 'o', 'Color', 'b', 'MarkerFaceColor', 'w', 'MarkerSize', 4);
plot(g_coh*100, g_pright_strg(:,1), 'b-');
plot(g_coh*100, g_pright_nostrg(:,1), 'b-.');
plot([g_coh(1) g_coh(end)]*100,[0.5 0.5], 'k--');
plot([0 0],[0 1],'k--');
set(gca, 'XLim', [-60 60], 'XTick', -60:10:60, 'XTickLabel', makeTickLabel(-60:10:60,20), ...
         'YLim',[-0.02 1.02], 'YTick', 0:0.1:1, 'YTickLabel', makeTickLabel(0:0.1:1,0.2), 'TickDir', 'out');
xlabel('Motion strength (%coh)');
ylabel('Probability right');
h = legend({'w/ stg','w/o strg'},'Location','NorthWest');
set(h,'FontSize',9);
legend('boxoff');
subplot(3,2,4); hold on;
errorbar(coh_axis/10, pright_strg(:,2), pright_strg_se(:,2), 'o', 'Color', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 4);
errorbar(coh_axis/10, pright_nostrg(:,2), pright_nostrg_se(:,2), 'o', 'Color', 'r', 'MarkerFaceColor', 'w', 'MarkerSize', 4);
plot(g_coh*100, g_pright_strg(:,2), 'r-');
plot(g_coh*100, g_pright_nostrg(:,2), 'r-.');
plot([g_coh(1) g_coh(end)]*100,[0.5 0.5], 'k--');
plot([0 0],[0 1],'k--');
set(gca, 'XLim', [-60 60], 'XTick', -60:10:60, 'XTickLabel', makeTickLabel(-60:10:60,20), ...
         'YLim',[-0.02 1.02], 'YTick', 0:0.1:1, 'YTickLabel', makeTickLabel(0:0.1:1,0.2), 'TickDir', 'out');
xlabel('Motion strength (%coh)');
ylabel('Probability right');

subplot(3,2,[5 6]); axis off;
filename2 = filename;
filename2(filename2=='_')= ' ';
if ~exist('modelLL','var')
    modelLL(2) = m2;
end

%         1  2 3     4      5    6  7    8        9      10        11     12     13
      %  [k  B theta stmeff bias dk dsig sig2beta weber guess_rate theta2 dtheta dtheta2] %%%%%%

text(0.1,0.3, sprintf( ... 
    'file: %s\n\t\t\t\t\t\tk= %f\n\t\t\t\t\t\tB= %f\n\t\t\t\t\t\ttheta= %f\n\t\t\t\t\t\tdcohstim= %f\n\t\t\t\t\t\tdcohbias= %f\n\t\t\t\t\t\tdk= %f\n\t\t\t\t\t\tdsigma2= %f\n\t\t\t\t\t\tsigma2beta= %f\n\t\t\t\t\t\tweber= %f\n\t\t\t\t\t\tguessrate= %f\n\t\t\t\t\t\ttheta2= %f\n\t\t\t\t\t\tdtheta= %f\n\t\t\t\t\t\tdtheta2= %f\n\nerr= %f\nerr2= %f\nnum trials = %d\ndur set = [%d %d]', ... 
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
    print('-dpdf',[folder filename '_dStoch_' fitwhat{1} '_' fitwhat{2} '_q=' num2str(q)]);
end

% ********** End, make the figure *********

end

if saveMAT
    save([folder filename '_dStoch_' fitwhat '.mat']);
end
