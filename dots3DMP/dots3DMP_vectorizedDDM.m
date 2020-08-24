
% SJ 06/23/2020
% vectorized version of DDM for sim_behavior..

% ...anyway, so far this does not provide a speedup on for loop
% seem to be some large overheads in bsxfun and randn

% I suppose for loop is fine to run over each trial once, but a vectorized
% solution for repeatedly simulating each trial (to generate likelihoods)
% would be good.
% presumably we also only need this once for each unique trial, and then an
% just lookup the likelihood for a particular trial configuration.


cohs = [0.2 0.5];
hdgs = [-10 -5 -2.5 -1.25 -eps eps 1.25 2.5 5 10];
mods = [1 2 3];
deltas = [-3 0 3];
nreps = 200;

duration = 2000;

kves = 1.2;
kvis = [0.8; 2];
B = 70;

sigmaVes = 1;
sigmaVis = [1 1];

Tnd = 300;

maxdur = duration;
k = mean([kves kvis']);

sigma = mean([sigmaVes sigmaVis]);
[~, ~, logOddsCorrMap, tAxis, vAxis] = makeLogOddsCorrMap_3DMP(hdgs,k,B,sigma,maxdur,0);

% create acceleration and velocity profiles (arbitrary for now)
% SJ 04/2020
% Hou et al. 2019, peak vel = 0.37m/s, SD = 210ms
vel = normpdf(1:duration,duration/2,210);
vel = 0.37*vel./max(vel);
acc = gradient(vel)*1000; % multiply by 1000 to get from m/s/ms to m/s/s

% normalize
vel = vel./max(vel);
acc = abs(acc./max(acc)); % (and abs)


%%

% vectorize
[hdg,modality,coh,delta,ntrials] = dots3DMP_create_trial_list(hdgs,mods,cohs,deltas,nreps);
dur = ones(ntrials,1) * duration;

tic
profile on
[~,cohI] = ismember(coh,cohs);
%cohsI    = sub2ind([length(coh),2],(1:length(coh))',cohI);

kvesV  = ones(size(hdg)).*kves;
kvisV  = kvis(cohI)'; 

% tic
% muVes = repmat(kvesV .* sind(hdg-delta/2),1,length(acc)) .* repmat(acc,length(kvesV),1);
% muVis = repmat(kvisV .* sind(hdg+delta/2),1,length(vel)) .* repmat(vel,length(kvisV),1);
% toc

muVes = bsxfun(@times,kvesV .* sind(hdg-delta/2), acc);
muVis = bsxfun(@times,kvisV .* sind(hdg+delta/2), vel);

wVes = sqrt(kvesV.^2 ./ (kvisV.^2 + kvesV.^2));
wVis = sqrt(kvisV.^2 ./ (kvisV.^2 + kvesV.^2));

profile off
profile viewer

%%
% mean drift rate
muComb = wVes.*muVes + wVis.*muVis;

muMat = nan(size(muVes));
muMat(modality==1,:) = muVes(modality==1,:);
muMat(modality==2,:) = muVis(modality==2,:);
muMat(modality==3,:) = muComb(modality==3,:);

% variance
sigmaVesV = ones(size(hdg)).*sigmaVes;
sigmaVisV = sigmaVis(cohI)';
sigmaComb = sqrt(wVes.^2 .* sigmaVesV.^2 + wVis.^2 .* sigmaVisV.^2); % assume zero covariance

sigmaMat = nan(size(hdg));
sigmaMat(modality==1) = sigmaVesV(modality==1);
sigmaMat(modality==2) = sigmaVisV(modality==2);
sigmaMat(modality==3) = sigmaComb(modality==3);

sigmaMat = repmat(sigmaMat,1,duration);

dv = [zeros(size(hdg)) cumsum(muMat + randn(size(muMat)).*sigmaMat,2)];

% dv = [zeros(size(hdg)) cumsum(normrnd(muMat,sigmaMat),2)];
% dvVes = [zeros(size(hdg)) cumsum(normrnd(muVes,sigmaVesV),2)];
% dvVis = [zeros(size(hdg)) cumsum(normrnd(muVis,sigmaVisV),2)];

% get RT, choice, and confidence from DV
[~,cRT] = max(abs(dv)>=B,[],2);
hitBound = cRT>1;

[RT,finalV] = deal(nan(size(cRT)));
RT(~hitBound) = dur(~hitBound) + Tnd; 
RT(hitBound)  = cRT(hitBound) + Tnd; 

durind = sub2ind(size(dv),(1:length(RT))',dur);
RTind  = sub2ind(size(dv),(1:length(RT))',cRT);

finalV(~hitBound) = dv(durind(~hitBound));
finalV(hitBound)  = B.*sign(dv(RTind(hitBound)));
choice = sign(finalV);

% use map to look up log-odds that the motion is rightward
[~,thisV] = min(abs(bsxfun(@minus,vAxis,finalV)),[],2);
[~,thisT] = min(abs(bsxfun(@minus,tAxis',RT)),[],2);

logmapinds = sub2ind(size(logOddsCorrMap),thisV,thisT);

logOddsCorr = logOddsCorrMap(logmapinds);
expectedPctCorr = logistic(logOddsCorr);
conf = 2*expectedPctCorr - 1;

choice(choice==0) = sign(randn); % not needed under usual circumstances
choice(choice==1) = 2; choice(choice==-1) = 1; % 1=left, 2=right


