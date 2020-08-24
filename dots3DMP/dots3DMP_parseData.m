%% parse data
% create and use matrices of summary data indexed by variables of interest

n = nan(length(mods),length(cohs),length(deltas)+1,length(hdgs));
                               % add extra column^ for pooling all trials irrespective of delta
pRight = n;
RTmean = n;
RTse = n;
confMean = n;
confSE = n;

xVals = hdgs(1):0.1:hdgs(end);
yVals = nan(length(mods),length(cohs),length(deltas)+1,length(xVals));

B = cell(length(mods),length(cohs),length(deltas)+1);
stats = cell(length(mods),length(cohs),length(deltas)+1);

for m = 1:length(mods)
for c = 1:length(cohs)
for d = 1:length(deltas)+1 % add extra column for all trials irrespective of delta

    for h = 1:length(hdgs)
        if d==length(deltas)+1
            J = data.modality==mods(m) & data.coherence==cohs(c) & data.heading==hdgs(h); % all trials irrespective of delta
        else
            J = data.modality==mods(m) & data.coherence==cohs(c) & data.heading==hdgs(h) & data.delta==deltas(d);
        end
        
        n(m,c,d,h) = sum(J);
        pRight(m,c,d,h) = sum(J & data.choice==2) / n(m,c,d,h); % 2 is rightward
        
        if RTtask
            RTmean(m,c,d,h) = mean(data.RT(J));
            RTse(m,c,d,h) = std(data.RT(J))/sqrt(n(m,c,d,h));
        end
        
        if conftask==1
            confMean(m,c,d,h) = mean(data.conf(J));
            confSE(m,c,d,h) = std(data.conf(J))/sqrt(n(m,c,d,h));
        else % PDW
            J = J & ~data.oneConfTargTrial;
            confMean(m,c,d,h) = sum(J & data.PDW==1) / sum(J); % 1 is high
            % SE gets calculated below
        end            
    end

    % fit logistic regression
    if d==length(deltas)+1
        K = data.modality==mods(m) & data.coherence==cohs(c); % all trials irrespective of delta
    else
        K = data.modality==mods(m) & data.coherence==cohs(c) & data.delta==deltas(d);
    end
    if sum(~isnan(data.heading(K)))>=3*length(hdgs) && length(hdgs)>5
        X = data.heading(K);
        y = data.choice(K)==2; % 2 is rightward
        [B{m,c,d}, ~, stats{m,c,d}] = glmfit(X, y, 'binomial');
        yVals(m,c,d,:) = glmval(B{m,c,d},xVals,'logit');
        plotLogistic(m,c,d) = 1;
    else
        plotLogistic(m,c,d) = 0;
    end

end
end
end

% standard error of proportion
pRightSE = sqrt( (pRight.*(1-pRight)) ./ n );
if conftask==2
    confSE = sqrt( (confMean.*(1-confMean)) ./ n );
end

% copy vestib-only data to both coherences, to aid plotting
n(1,2,:,:) = n(1,1,:,:);
pRight(1,2,:,:) = pRight(1,1,:,:);
pRightSE(1,2,:,:) = pRightSE(1,1,:,:);
confMean(1,2,:,:) = confMean(1,1,:,:);
confSE(1,2,:,:) = confSE(1,1,:,:);
RTmean(1,2,:,:) = RTmean(1,1,:,:);
RTse(1,2,:,:) = RTse(1,1,:,:);
confMean(1,2,:,:) = confMean(1,1,:,:);
confSE(1,2,:,:) = confSE(1,1,:,:);
yVals(1,2,:,:) = yVals(1,1,:,:);
plotLogistic(1,2,:) = plotLogistic(1,1,:);
B(1,2,:) = B(1,1,:);
stats(1,2,:) = stats(1,1,:);

% some datasets have three cohs!
if length(cohs)==3
    n(1,3,:,:) = n(1,1,:,:);
    pRight(1,3,:,:) = pRight(1,1,:,:);
    pRightSE(1,3,:,:) = pRightSE(1,1,:,:);
    confMean(1,3,:,:) = confMean(1,1,:,:);
    confSE(1,3,:,:) = confSE(1,1,:,:);
    RTmean(1,3,:,:) = RTmean(1,1,:,:);
    RTse(1,3,:,:) = RTse(1,1,:,:);
    confMean(1,3,:,:) = confMean(1,1,:,:);
    confSE(1,3,:,:) = confSE(1,1,:,:);
    yVals(1,3,:,:) = yVals(1,1,:,:);
    plotLogistic(1,3,:) = plotLogistic(1,1,:);
    B(1,3,:) = B(1,1,:);
    stats(1,3,:) = stats(1,1,:);
end

if length(cohs)>3
    error('can''t handle >3 cohs');
end

