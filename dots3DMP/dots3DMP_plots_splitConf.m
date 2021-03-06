% plot curves split by high/low confidence
% bit of a kluge, just splitting by median confidence rating, but eventual
% plotting will depend on exact analyses, so just use this as a sanity
% check that confidence ratings are genuinely saying something about
% perceived p(correct)

% first, for all trials irrespective of delta
D = length(deltas)+1; % (the extra column we made for pooling across deltas)
% OR select just delta=0:
% D = find(deltas==0);

         %ves %vis %comb
clr{1} = {'ko','mo','co'};
clr{2} = {'ko','ro','bo'};
clr{3} = {'ko','yo','go'};
figure(111+D);
set(gcf,'Color',[1 1 1],'Position',[300 500 950+300*(length(cohs)-2) 800],'PaperPositionMode','auto'); clf;
for c = 1:length(cohs)
    subplot(2+(~isnan(RTmean(1,1,2))),length(cohs),c);
    for m = 1:length(mods)     % m c d h
        if plotLogistic(m,c,D)
            h1(m) = plot(xVals,squeeze(yVals(m,c,D,:,1)),[clr{c}{m}(1) '-']); hold on;
            h2(m) = plot(xVals,squeeze(yVals(m,c,D,:,2)),[clr{c}{m}(1) '--']); hold on;

            %errorbar(hdgs, squeeze(pRight(m,c,D,:)), squeeze(pRightSE(m,c,D,:)), clr{c}{m});
        else
            %h(m) = errorbar(hdgs, squeeze(pRight(m,c,D,:)), squeeze(pRightSE(m,c,D,:)), [clr{c}{m} '-']); hold on;
        end
        ylim([0 1]);
        if length(mods)>1; title(['coh = ' num2str(cohs(c))]); end
    end
    legend(h1,'vestib','visual','comb','Location','northwest');
    xlabel('heading angle (deg)'); ylabel('proportion rightward choices');

    subplot(2+(~isnan(RTmean(1,1,2))),length(cohs),c+length(cohs));
    for m = 1:length(mods)
        h1(m) = errorbar(hdgs, squeeze(confMean(m,c,D,:,1)), squeeze(confSE(m,c,D,:,1)), [clr{c}{m} '-']);
        h2(m) = errorbar(hdgs, squeeze(confMean(m,c,D,:,2)), squeeze(confSE(m,c,D,:,2)), [clr{c}{m} '--']);

        ylim([0 1.6]); hold on;
    end
    xlabel('heading angle (deg)'); ylabel('saccadic endpoint (''confidence'', %)');
    
    if ~isnan(RTmean(1,1,2))
        subplot(3,length(cohs),c+length(cohs)*2);
        for m = 1:length(mods)
            h1(m) = errorbar(hdgs, squeeze(RTmean(m,c,D,:,1)), squeeze(RTse(m,c,D,:,1)), [clr{c}{m} '-']); hold on;
            h2(m) = errorbar(hdgs, squeeze(RTmean(m,c,D,:,2)), squeeze(RTse(m,c,D,:,2)), [clr{c}{m} '--']); hold on;

        end
        xlabel('heading angle (deg)'); ylabel('RT (s)');
    end
    
end