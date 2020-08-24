%% parse data
% create and use matrices of summary data indexed by variables of interest

n = nan(length(mods),length(cohs),length(deltas)+1,length(hdgs));
                               % add extra column^ for pooling all trials irrespective of delta
[pRight, RTmean, RTse, confMean, confSE] = deal(n);

xVals = hdgs(1):0.1:hdgs(end);
yVals = nan(length(mods),length(cohs),length(deltas)+1,length(xVals));

B = cell(length(mods),length(cohs),length(deltas)+1);
stats = cell(length(mods),length(cohs),length(deltas)+1);

% split confidence into high and low for saccadic endpoint data 
% let's start just by median...
% eventually, probably need to split each subject separately (or normalise
% within each subject to have a single, equivalent, scale)
% easier for PDW

if conftask==1
    hiConf = data.conf >= median(data.conf);
elseif conftask==2 
    hiConf = data.PDW;
end

for m = 1:length(mods)
for c = 1:length(cohs)
for d = 1:length(deltas)+1 % add extra column for all trials irrespective of delta

    for h = 1:length(hdgs)
        if d==length(deltas)+1
            J = data.modality==mods(m) & data.coherence==cohs(c) & data.heading==hdgs(h); % all trials irrespective of delta
        else
            J = data.modality==mods(m) & data.coherence==cohs(c) & data.heading==hdgs(h) & data.delta==deltas(d);
        end
        Jhi = J & hiConf;
        Jlo = J & ~hiConf;
        
        n(m,c,d,h,1) = sum(Jhi);
        n(m,c,d,h,2) = sum(Jlo);

        pRight(m,c,d,h,1) = sum(Jhi & data.choice==2) / n(m,c,d,h,1); % 2 is rightward
        pRight(m,c,d,h,2) = sum(Jlo & data.choice==2) / n(m,c,d,h,2); % 2 is rightward

        if RTtask
            RTmean(m,c,d,h,1) = mean(data.RT(Jhi));
            RTmean(m,c,d,h,2) = mean(data.RT(Jlo));
            
            RTse(m,c,d,h,1) = std(data.RT(Jhi))/sqrt(n(m,c,d,h,1));
            RTse(m,c,d,h,2) = std(data.RT(Jlo))/sqrt(n(m,c,d,h,2));
        end
        
        if conftask==1
            confMean(m,c,d,h,1) = mean(data.conf(Jhi));
            confSE(m,c,d,h,1) = std(data.conf(Jhi))/sqrt(n(m,c,d,h,1));
            
            confMean(m,c,d,h,2) = mean(data.conf(Jlo));
            confSE(m,c,d,h,2) = std(data.conf(Jlo))/sqrt(n(m,c,d,h,2));
            
        else % PDW
            confMean(m,c,d,h,1) = sum(J & data.PDW==1) / n(m,c,d,h,1); % 1 is high
            confMean(m,c,d,h,2) = sum(J & data.PDW==0) / n(m,c,d,h,2); 

            % SE gets calculated below
        end            
    end

    % fit logistic regression
    if d==length(deltas)+1
        K = data.modality==mods(m) & data.coherence==cohs(c); % all trials irrespective of delta
    else
        K = data.modality==mods(m) & data.coherence==cohs(c) & data.delta==deltas(d);
    end
    Khi = K & hiConf;

    if sum(~isnan(data.heading(Khi)))>=3*length(hdgs) && length(hdgs)>5
        X = data.heading(Khi);
        y = data.choice(Khi)==2; % 2 is rightward
        [B{m,c,d,1}, ~, stats{m,c,d,1}] = glmfit(X, y, 'binomial');
        yVals(m,c,d,:,1) = glmval(B{m,c,d,1},xVals,'logit');
        plotLogistic(m,c,d,1) = 1;
    else
        plotLogistic(m,c,d,1) = 0;
    end
    
    Klo = K & ~hiConf;
    if sum(~isnan(data.heading(Klo)))>=3*length(hdgs) && length(hdgs)>5
        X = data.heading(Klo);
        y = data.choice(Klo)==2; % 2 is rightward
        [B{m,c,d,2}, ~, stats{m,c,d,2}] = glmfit(X, y, 'binomial');
        yVals(m,c,d,:,2) = glmval(B{m,c,d,2},xVals,'logit');
        plotLogistic(m,c,d,2) = 1;
    else
        plotLogistic(m,c,d,2) = 0;
    end

end
end
end

pRightSE = sqrt( (pRight.*(1-pRight)) ./ n );
if conftask==2
    confSE = sqrt( (confMean.*(1-confMean)) ./ n );
end

% copy vestib-only data to both coherences, to aid plotting
n(1,2,:,:,:) = n(1,1,:,:,:);
pRight(1,2,:,:,:) = pRight(1,1,:,:,:);
pRightSE(1,2,:,:,:) = pRightSE(1,1,:,:,:);
confMean(1,2,:,:,:) = confMean(1,1,:,:,:);
confSE(1,2,:,:,:) = confSE(1,1,:,:,:);
RTmean(1,2,:,:,:) = RTmean(1,1,:,:,:);
RTse(1,2,:,:,:) = RTse(1,1,:,:,:);
confMean(1,2,:,:,:) = confMean(1,1,:,:,:);
confSE(1,2,:,:,:) = confSE(1,1,:,:,:);
yVals(1,2,:,:,:) = yVals(1,1,:,:,:);
plotLogistic(1,2,:,:) = plotLogistic(1,1,:,:);
B(1,2,:,:) = B(1,1,:,:);
stats(1,2,:,:) = stats(1,1,:,:);

% some datasets have three cohs!
if length(cohs)==3
    n(1,3,:,:,:) = n(1,1,:,:,:);
    pRight(1,3,:,:,:) = pRight(1,1,:,:,:);
    pRightSE(1,3,:,:,:) = pRightSE(1,1,:,:,:);
    confMean(1,3,:,:,:) = confMean(1,1,:,:,:);
    confSE(1,3,:,:,:) = confSE(1,1,:,:,:);
    RTmean(1,3,:,:,:) = RTmean(1,1,:,:,:);
    RTse(1,3,:,:,:) = RTse(1,1,:,:,:);
    confMean(1,3,:,:,:) = confMean(1,1,:,:,:);
    confSE(1,3,:,:,:) = confSE(1,1,:,:,:);
    yVals(1,3,:,:,:) = yVals(1,1,:,:,:);
    plotLogistic(1,3,:,:) = plotLogistic(1,1,:,:);
    B(1,3,:,:) = B(1,1,:,:);
    stats(1,3,:,:) = stats(1,1,:,:);
end

if length(cohs)>3
    error('can''t handle >3 cohs');
end

