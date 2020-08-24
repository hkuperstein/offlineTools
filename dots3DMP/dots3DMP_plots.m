%% first, for all trials irrespective of delta
D = length(deltas)+1; % (the extra column we made for pooling across deltas)
% OR select just delta=0:
% D = find(deltas==0);

         %ves %vis %comb
clr{1} = {'ko','mo','co'};
clr{2} = {'ko','ro','bo'};
clr{3} = {'ko','yo','go'};
figure(101+D);
set(gcf,'Color',[1 1 1],'Position',[300 500 950+300*(length(cohs)-2) 800],'PaperPositionMode','auto'); clf;
for c = 1:length(cohs)
    subplot(2+(~isnan(RTmean(1,1,2))),length(cohs),c);
    for m = 1:length(mods)     % m c d h
        if plotLogistic(m,c,D)
            h(m) = plot(xVals,squeeze(yVals(m,c,D,:)),[clr{c}{m}(1) '-']); hold on;
            errorbar(hdgs, squeeze(pRight(m,c,D,:)), squeeze(pRightSE(m,c,D,:)), clr{c}{m});
        else
            h(m) = errorbar(hdgs, squeeze(pRight(m,c,D,:)), squeeze(pRightSE(m,c,D,:)), [clr{c}{m} '-']); hold on;
        end
        ylim([0 1]);
        if length(mods)>1; title(['coh = ' num2str(cohs(c))]); end
    end
    legend(h,'vestib','visual','comb','Location','northwest');
    xlabel('heading angle (deg)'); ylabel('proportion rightward choices');

    subplot(2+(~isnan(RTmean(1,1,2))),length(cohs),c+length(cohs));
    for m = 1:length(mods)
        h(m) = errorbar(hdgs, squeeze(confMean(m,c,D,:)), squeeze(confSE(m,c,D,:)), [clr{c}{m} '-']);
        ylim([0 1]); hold on;
    end
    xlabel('heading angle (deg)'); ylabel('saccadic endpoint (''confidence'', %)');
    
    if ~isnan(RTmean(1,1,2))
        subplot(3,length(cohs),c+length(cohs)*2);
        for m = 1:length(mods)
            h(m) = errorbar(hdgs, squeeze(RTmean(m,c,D,:)), squeeze(RTse(m,c,D,:)), [clr{c}{m} '-']); hold on;
        end
        xlabel('heading angle (deg)'); ylabel('RT (s)');
    end
    
end


%% now separate by delta

if length(deltas)>1

clr{1} = {'bs','cs','gs'};
clr{2} = {'b^','c^','g^'};
clr{3} = {'bo','co','go'};

clear L;
figure(108);
set(gcf,'Color',[1 1 1],'Position',[50 20 950+300*(length(cohs)-2) 800],'PaperPositionMode','auto'); clf;
for c = 1:length(cohs)
    subplot(2+(~isnan(RTmean(1,1,2))),length(cohs),c);
    for d = 1:length(deltas)     % m c d h
        if plotLogistic(3,c,d)
            h(d) = plot(xVals,squeeze(yVals(3,c,d,:)),[clr{c}{d}(1) '-']); hold on;
            errorbar(hdgs, squeeze(pRight(3,c,d,:)), squeeze(pRightSE(3,c,d,:)), clr{c}{d});
        else
            h(d) = errorbar(hdgs, squeeze(pRight(3,c,d,:)), squeeze(pRightSE(3,c,d,:)), [clr{c}{d} '-']); hold on
        end
        L{d} = sprintf('\x0394=%d',deltas(d));
        ylim([0 1]);
        if length(mods)>1; title(['coh = ' num2str(cohs(c))]); end
    end
    legend(h,L,'location','northwest');
    xlabel('heading angle (deg)'); ylabel('proportion rightward choices');

    subplot(2+(~isnan(RTmean(1,1,2))),length(cohs),c+length(cohs));
    for d = 1:length(deltas)
        h(d) = errorbar(hdgs, squeeze(confMean(3,c,d,:)), squeeze(confSE(3,c,d,:)), [clr{c}{d} '-']);
        ylim([0 1]); hold on;
    end
    xlabel('heading angle (deg)'); ylabel('saccadic endpoint (''confidence'', %)');
    
    if ~isnan(RTmean(1,1,2))
        subplot(3,length(cohs),c+length(cohs)*2);
        for d = 1:length(deltas)
            h(d) = errorbar(hdgs, squeeze(RTmean(3,c,d,:)), squeeze(RTse(3,c,d,:)), [clr{c}{d} '-']); hold on;
        end
        xlabel('heading angle (deg)'); ylabel('RT (s)');
    end
    
end

end