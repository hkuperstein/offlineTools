function [X, err_final, fit, fitInterp] = dots3DMP_fitDDM(data,options,guess,fixed)


if options.ploterr
    options.fh = figure;
    hold on;   
    xlabel('Call number' );
    ylabel('err');
end


% parameter bounds for fitting
%      [kves kvisMult B]
% X_lo = [0.4 1 20 ];
    X_lo = guess/3;
% X_hi = [3   9 120];
    X_hi = guess*3;



global call_num; call_num=1;

if all(fixed)
    X = guess;
else

switch options.fitMethod
    case 'fms'
        fitOptions = optimset('Display', 'final', 'MaxFunEvals', 500*sum(fixed==0), 'MaxIter', ... 
            500*sum(fixed==0), 'TolX', 1e-3, 'TolFun', 1e-2, 'UseParallel', 'Always');
        [X, ~] = fminsearch(@(x) dots3DMP_fitDDM_err(x,data,options), guess, fitOptions);

    case 'global'
        % GlobalSearch from Global Optimization Toolbox
        fitOptions = optimoptions(@fmincon,'Display','iter',...
            'Algorithm','interior-point',...
            'FinDiffType','central',...
            'UseParallel','always');
        problem = createOptimProblem('fmincon','x0',guess(fixed==0),'objective',...
            @(x) dots3DMP_fitDDM_err(x,data,options),'lb',X_lo(fixed==0),...
            'ub',X_hi(fixed==0),'options',fitOptions);
        gs = GlobalSearch;
        [X,~,~,~,~] = run(gs,problem);      

    case 'multi'
        % MultiStart from Global Optimization Toolbox
        fitOptions = optimoptions(@fmincon,'Display','iter',...
            'Algorithm','interior-point',...
            'FinDiffType','central',...
            'UseParallel','always');
        problem = createOptimProblem('fmincon','x0',guess(fixed==0),'objective',...
            @(x) dots3DMP_fitDDM_err(x,data,options),'lb',X_lo(fixed==0),...
            'ub',X_hi(fixed==0),'options',fitOptions);
        ms = MultiStart('UseParallel','always');
        [X,~,~,~,~] = run(ms,problem,200);    
    
    case 'pattern'
        % Pattern Search from Global Optimization Toolbox
        fitOptions = psoptimset(@patternsearch);
        fitOptions = psoptimset(fitOptions,'CompletePoll','on','Vectorized','off',...
            'UseParallel','always');
        [X,~,~,~] = patternsearch(@(x) dots3DMP_fitDDM_err(x,data,options),...
            guess(fixed==0),[],[],[],[],X_lo(fixed==0),X_hi(fixed==0),[],fitOptions);
end

end

% run err func again at the fitted/fixed params to generate a final
% error value and model-generated data points (trial outcomes)
options.ploterr = 0;
[err_final, fit] = dots3DMP_fitDDM_err(X,data);

% *** now run it one more time with interpolated headings
% to generate smooth plots ***
    % (this ordinarily works well to show the fits vs. data, but in the
    % current version the model fits/predictions are generated by monte
    % carlo simulation and are thus too noisy to be worth repeating at
    % finely spaced headings (you can try uncommenting the next two lines
    % to see what this looks like). So for now hdgs will be the actual hdgs\
    % from data and this section will seem wholly redundant, but will be
    % useful later)

% nsteps = 33; % should be odd so there's a zero
% hdgs = linspace(min(data.heading),max(data.heading),nsteps);
hdgs = unique(data.heading)'; % TEMP, see above comment

% Dfit is a dummy dataset with the same proportions of all trial types,
% just repeated for the new interpolated headings (and some arbitrary nreps)
mods = unique(data.modality);
cohs = unique(data.coherence);
deltas = unique(data.delta);

% /begin generate condition list (could be moved into a function; it's the
% same as the sim code and maybe elsewhere)
nreps = 200;

% repeat heading list once for ves, ncoh for vis, and ncoh x ndelta for comb
numHdgGroups = any(ismember(mods,1)) + ...
               any(ismember(mods,2)) * length(cohs) + ...
               any(ismember(mods,3)) * length(cohs)*length(deltas);
hdg = repmat(hdgs', numHdgGroups, 1);

lh = length(hdgs);
coh = nan(size(hdg));
delta = nan(size(hdg));
modality = nan(size(hdg));

% kluge for ves, call it the lowest coh by convention
if any(ismember(mods,1))
    coh(1:lh) = cohs(1); 
    delta(1:lh) = 0;
    modality(1:lh) = 1;
    last = lh;
else
    last = 0;    
end

if any(ismember(mods,2)) % loop over coh for vis
    for c = 1:length(cohs)
        these = last+(c-1)*lh+1 : last+(c-1)*lh+lh;
        coh(these) = cohs(c);
        delta(these) = 0;
        modality(these) = 2;
    end
    last = these(end);
end

if any(ismember(mods,3)) % loop over coh and delta for comb
    for c = 1:length(cohs)
        for d = 1:length(deltas)
            here = last + 1 + (c-1)*lh*length(deltas) + (d-1)*lh;
            these = here:here+lh-1;
            coh(these) = cohs(c);
            delta(these) = deltas(d);
            modality(these) = 3;
        end
    end
end

% now replicate times nreps
condlist = [hdg modality coh delta];
trialTable = repmat(condlist,nreps,1);
Dfit.heading = trialTable(:,1);  
Dfit.modality = trialTable(:,2);  
Dfit.coherence = trialTable(:,3);  
Dfit.delta = trialTable(:,4);

%/end generate condition list

% the observables are just placeholders because the err func still needs to
% calculate err even though in this case it's meaningless. What will be
% plotted is fitInterp, not Dfit
Dfit.choice = ones(size(Dfit.heading));
Dfit.RT = ones(size(Dfit.heading));
Dfit.conf = ones(size(Dfit.heading));

[~,fitInterp] = dots3DMP_fitDDM_err(X,Dfit);


end



