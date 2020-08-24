% offline analysis wrapper for PLDAPS data, Dots protocols
% CF spawned it from dots3DMP_offlineAnalysis, begun 10-3-18

clear all; close all

% these will specify which (previously saved) mat file to load
subject = 'hanzo';
dateRange = 20190401:20200315; % must match previously saved mat file

folder = '/Users/chris/Documents/MATLAB/PLDAPS_data/';
file = [subject '_' num2str(dateRange(1)) '-' num2str(dateRange(end)) '.mat'];
load([folder file], 'data');

% remove invalid trials (fixation breaks, etc)
removethese = isnan(data.choice) | isnan(data.PDW); % may need to add conditions here
fnames = fieldnames(data);
for F = 1:length(fnames)
    eval(['data.' fnames{F} '(removethese) = [];']);
end


% define signed coherence; for now positive = rightward (0 deg)
data.scoh = data.coherence;
data.scoh(data.direction==180) = -data.scoh(data.direction==180);


%% parse data

cohs = unique(data.scoh);
dirs = unique(data.direction);

n = nan(length(cohs),1);
pRight = n;
pHigh = n;
RTmean = n;
RTse = n;

% for logistic regression (or whatever)
xVals = linspace(cohs(1),cohs(end),100);

for c = 1:length(cohs)
    J = data.scoh==cohs(c);
    n(c) = sum(J);
    pRight(c) = sum(J & data.choice==2) / n(c); % 1 is left, 2 is right
    pHigh(c) = sum(J & data.PDW==1) / n(c); % 1 is high-bet
    RTmean(c) = mean(data.RT(J));
    RTse(c) = std(data.RT(J))/sqrt(n(c));
end
pRightSE = sqrt( (pRight.*(1-pRight)) ./ n );
pHighSE = sqrt( (pHigh.*(1-pHigh)) ./ n );

% fit logistic regression
X = data.scoh;
y = data.choice==2; % 2 is right
[B, ~, stats] = glmfit(X, y, 'binomial');
yVals = glmval(B,xVals,'logit');


%% plot

figure(101);
set(gcf,'Color',[1 1 1],'Position',[300 500 450 800],'PaperPositionMode','auto'); clf;
subplot(2,1,1);
plot(xVals,yVals,'b-'); hold on;
errorbar(cohs, pRight, pRightSE, 'bo');
set(gca,'xtick',cohs([1 2 4 6 8 10 11]),'tickdir','out');
ylim([0 1]); xlim([-0.55 0.55]);
xlabel('motion strength (% coh)'); ylabel('proportion rightward choices');
changeAxesFontSize(gca,14,14);

subplot(2,1,2);
errorbar(cohs, pHigh, pHighSE, 'ro');
set(gca,'xtick',cohs([1 2 4 6 8 10 11]),'tickdir','out');
ylim([0 1]); xlim([-0.55 0.55]);
xlabel('motion strength (% coh)'); ylabel('proportion high bet');
changeAxesFontSize(gca,14,14);








