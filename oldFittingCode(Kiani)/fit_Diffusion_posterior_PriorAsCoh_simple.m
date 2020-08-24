function [modelParam,modelLL,trialData] = ...
            fit_Diffusion_posterior_PriorAsCoh(trialData, guess, fixed, fitwhat, orig_coh_set, orig_coh_freq, filename)
%
% function [modelParam,modelLL,trialData] = 
%               fit_Diffusion_posterior_PriorAsCoh(trialData, guess, fixed, fitwhat, orig_coh_set, orig_coh_freq)
% 
% fits probability of rightward and sure target choices on stimulated and
% non-stimulated trials using a diffusion model. 
% 
% The model has the following free parameters
%   K               the coefficient that converts the motion strength to evidence.
%   B               bound height.
%   logodds_crit    the criterion on log-odds correct that governs choosing the sure target. 
%   weber_fraction  potentially useful, but turns out to be ineffective here. 
%   uStim_eff       the size of the micto-stim effect expressed in effective coherence.
%   init_bias       the size of bias in non-stimulated trials. for fitting purposes it
%                   may be better to get the value from the data and keep this parameter
%                   fixed.
% 
% Input:
% 
% trialData     a structure made from FIRA, see uStim-Fit.m
% 
% guess     initial guess for the model parameters
% 
% fixed     a vector that defines which of the model parameters should be treated as fixed parameters 
%   
% fitwhat   is a cell with two elements and defines what should be fitted for
%           non-stimulated and stimulated trials. the two elements govern fitting of
%           non-stimulated and stimulated trials, respectively.
%           each element can be 
%                  'PrightPS'   fit both probability right and probability sure target. 
%                  'Pright'     fit only probability right. uses trials in which sure 
%                           target was not shwon and the trials in which sure target
%                           was shown but waived.
%                  'PS'     fit only probability sure target. fits Pright for trials
%                           without sure target and Pstrg for trials in which sure
%                           target was shown.
%                   'Pright-woTs' fit only Pright and only for trials in which Ts was
%                           not shown.
%           playing with the two elements allows a variety of fitting strategies. for
%           example you can choose {'PrightPS','Pright'} which fits all the free
%           parameters based on all non-stimulated trials and the stimulated trials
%           in which sure target was waived. Then you can predict probability of
%           choosing sure target.
% 
% orig_coh_set      defines the set of coherences that should be used for the calculation
%                   of marginals. leave it empty if all the coherences in trialData
%                   should be used.
% 
% orig_coh_freq     prior probability of cohereces defined in orig_coh_set
% 
% 
% RK, 10/2009
% 

    dbstop if error
    expectedPright = [];
    expectedPS = [];
%     global callNum_;
    callNum_ = 1;
    global err2;
    
    modelParam = struct('init' ,[] ,'fixed', [] ,'final', []);
    
    if nargin<6
        orig_coh_set = [];
        orig_coh_freq = [];
    end
    
    if nargin < 5 || isempty(fitwhat)
        fitwhat = {'PrightPS', 'PrightPS'};
    end
    
%         get the coherence levels, the frequency of each coherence level and the maximum
%         stimulus duration 
%         Note that coh must be signed for the algorithm to work properly 
    coh_set = unique(trialData.Coh);
    coh_set_freq = nan(size(coh_set));
    for cc = 1 : length(coh_set)
        coh_set_freq(cc) = sum(trialData.Coh==coh_set(cc))/length(trialData.Coh);
    end
    max_dur = max(trialData.StimDur);

        %do the fitting
    modelParam.init = guess;
    modelParam.fixed = fixed;
    modelParam.final = guess;
    if all(fixed==1)
        fprintf ( '\n\nall parameters are fixed\n\n' );
%         [expectedPright, expectedPS, err] = fitmodel_MLEerr_done(modelParam.final, guess, fixed, trialData, fitwhat, orig_coh_set, orig_coh_freq, filename);
        err = fitmodel_MLEerr(modelParam.final);
        modelLL = -err;
        trialData.expectedPright = expectedPright;
        trialData.expectedPS = expectedPS;
        return;
    end
    
    options = optimset('Display' ,'final' ,'MaxFunEvals' ,500*sum(fixed==0) ,'MaxIter' ,500*sum(fixed==0), 'TolX', 1e-3, 'TolFun', 1e-2);
%    options = optimset('Display', 'final', 'MaxFunEvals', 1, 'MaxIter', 1, 'TolX', 0.1, 'TolFun', 0.1);
   
%     % termination tolerance on x is measured against the diameter of an n+1-dimensional
%     % simplex, meaning that the params should be roughly the same order of magnitude or else the
%     % larger ones will dominate the minimization.  Thus, normalize some of the params as follows...
%     guess(2) = guess(2)/100;
%     guess(3) = guess(3)/2;
%     guess(6) = guess(6)*10;
            % abandoned: did not seem to matter -CF
   
    [p,fval,exitflag] = fminsearch(@(x) fitmodel_MLEerr(x), guess(fixed==0), options);

    % sometimes worth trying fminunc instead
%     [p,fval,exitflag,~,~,hessian] = fminunc(@(x) fitmodel_MLEerr(x), guess(fixed==0), options);
%     modelParam.hessian = hessian;
%     modelParam.SEs = sqrt(diag(inv(hessian)));
    
    modelParam.final = getParam(p,guess,fixed);
    modelParam.exitflag = exitflag;

%     % ...and then re-scale them back to normal
%     modelParam.final(2) = modelParam.final(2)*100;
%     modelParam.final(3) = modelParam.final(3)*2;
%     modelParam.final(6) = modelParam.final(6)/10;
                % abandoned: did not seem to matter -CF

    modelLL = -fval;
    
        %calculated expectedPright and expectedPS for the fit parameters 
    guess = modelParam.final;
    fixed = ones(size(guess));
    fitmodel_MLEerr([]);    
    %[expectedPright, expectedPS, err] = fitmodel_MLEerr_done(modelParam.final, guess, fixed, trialData, fitwhat, orig_coh_set, orig_coh_freq, filename);
    trialData.expectedPright = expectedPright;
    trialData.expectedPS = expectedPS;
    
    % end of fitting
    
    
    function err = fitmodel_MLEerr(param)
        
        
        boundtype = 'flat';
%         boundtype = 'hyperbolic';       
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The model has the following free parameters
        %   k   the factor that converts the motion strength to evidence   e=kC
        %   B   the bound height
        %   logodds_crit ("theta")
        %   weber_fraction (not used)
        %   uStim_eff
        %   init_bias
        %   thetaOffset (rarely used)
        %   stim_choice_bias (old, never used)
        %   stim_k_mult
        %   stim_B_mult
        %   varOffset (conditional on uStim, like k_mult and B_mult)
        %   coh_var_term ("Beta")
        %   coh_var_offset (uStim offset to beta; not used)
        
        
        init_bias_influences_marginals = 0;
        stim_bias_influences_marginals = 0;

        param = getParam(param);

        % % param scaling (see calling function)
        % param(2) = param(2)*100;
        % param(3) = param(3)*2;
        % param(6) = param(6)/10;
                    % abandoned: did not seem to matter -CF

        k = param(1);        
        B = 30; % at some point fix it to 1
%         B = abs(param(2));  %don't accept negative bound heights
        sigma = abs(param(2));  %don't accept negative sigmas

        logodds_crit = abs(param(3));
        weber_fraction = abs(param(4));
        uStim_eff = param(5);
        init_bias = param(6);
        thetaOffset = param(7);
        if numel(param)>7
            stim_choice_bias = param(8);
        else
            stim_choice_bias = 0;
        end
        if numel(param)>8
            stim_k_mult = param(9);
        else
            stim_k_mult = 1;
        end
        if numel(param)>9
            stim_B_mult = param(10);
        else
            stim_B_mult = 1;
        end
        if numel(param)>10
            varOffset = param(11);
        else
            varOffset = 0;
        end
            % this is ugly, but ignore bound params and replace 12+13 with new variance params
%         if numel(param)>11
%             uInf = param(12);
%         else
            uInf = NaN;
            boundtype = 'flat';
%         end
%         if numel(param)>12
%             tau05 = param(13);
%         else
            tau05 = NaN;
%             boundtype = 'flat';
%         end
        if numel(param)>11
            coh_var_term = param(12);
        else
            error('cannot go back to old parameterization by omitting param 12');
        end
        if numel(param)>12
            coh_var_offset = param(13);
        else
            coh_var_offset = 0;
        end
        % done w/ setting free params
        
            %define the grid resolution for time and decision variable
        dt = 0.1;
        dx = min(0.1, B/100);
                
            %define the time axis 
        t = dt:dt:max_dur;
        
            %define xmesh and the ghost cells for bound crossing probability
        b_margin = [4*sigma; 4*sigma];
        xmesh = (-B-b_margin(1)+dx : dx : B+b_margin(2)-dx)';
            % adjust the mesh so that it includes zero and is symmetric 
        while ~any(xmesh==0)
            [~, I] = min(abs(xmesh));
            delta_mesh = xmesh(I);
            b_margin = b_margin + delta_mesh;
            xmesh = (-B-b_margin(1)+dx : dx : B+b_margin(2)-dx)'; % recompute the mesh
        end
        if mod(length(xmesh),2)~=1
            error('the length of xmesh must be odd');
        end

            %define a delta function on xmesh
        delta = zeros(size(xmesh));
        delta(abs(xmesh)==min(abs(xmesh))) = 1;

            %find the probability density of decision variable across time for stimulated and
            %nonstimulated trials 
        Ptb_coh = zeros(max_dur, 2, 2, length(coh_set));                % time * bound * stim * signed_coh
        Pxt_coh = zeros(length(xmesh), max_dur, 2, length(coh_set));    % xmesh * time * stim * signed_coh

        for s = 0 : 1
                %do not run FP4 unless it is really necessary 
            if stim_bias_influences_marginals==0
                if init_bias_influences_marginals==0 && s==0 && all(trialData.uStim==1)
                    continue;
                elseif s==1 && all(trialData.uStim==0)
                    continue;
                end
            end
                %run FP4 to calculate Ptb and Pxt for each coherence and stimulation condition 
            for c = 1 : length(coh_set)
                if s==1 && uStim_eff==0     %do not need to recalculate if uStim_eff is zero
                    Ptb_coh(:,:,s+1,c) = Ptb_coh(:,:,1,c);
                    Pxt_coh(:,:,s+1,c) = Pxt_coh(:,:,1,c);
                else
                    % must re-initialize b_change, which changes after each call of FP4
                    [b_change, ~] = bcinit(B,t,dt,boundtype,uInf,tau05);
                    if s==1
                        mu = k * stim_k_mult * (coh_set(c)+s*uStim_eff) + init_bias;
                        b_change(:,1) = b_change(:,1)*min(stim_B_mult,1) - (min(stim_B_mult,1)*B-B); % implements stim_B_mult using
                        b_change(:,2) = b_change(:,2)*min(stim_B_mult,1) + (min(stim_B_mult,1)*B-B); % existing b_change functionality
                    else
                        mu = k * (coh_set(c)+s*uStim_eff) + init_bias;
                    end
                    uinit = delta;

                    % new expression for variance
                    variance = (sigma^2+varOffset*s) + coh_var_term*abs(coh_set(c));
                    
                      [~, ~, Ptb_, ~, Pxt_] = FP4(xmesh, uinit, mu, sqrt(variance), b_change, b_margin, dt);

                        %store the arrays for each coherence level, use 1 ms time resolution
                    Ptb_coh(:,1,s+1,c) = local_sum(Ptb_(2:end,1),1/dt);  %lower bound
                    Ptb_coh(:,2,s+1,c) = local_sum(Ptb_(2:end,2),1/dt);  %upper bound
                    Pxt_coh(:,:,s+1,c) = Pxt_(:,1/dt:1/dt:end);
                end
            end
        end
        
            %find the region in which sure target should be chosen, use only
            %non-stimulated trials
        Ptb_marginal = zeros(max_dur, 2, 2);                            % time * bound * motion_direction(1 right, 2 left)
        Pxt_marginal = zeros(length(xmesh), max_dur, 2);                % xmesh * time * motion_direction

            %marginalize Ptb and Pxt across coherence levels for rightward and leftward motion 
        for c = 1 : length(coh_set)

            if ~isempty(orig_coh_set)
                [e,M] = ismember(coh_set(c), orig_coh_set);
                if e==1
                    F = orig_coh_freq(M);
                else
                    continue;
                end
            else
                F = coh_set_freq(c);
            end

                %divide F by 2 if both non-stimulated and stimulated trials are used to
                %calculate the marginals 
            if stim_bias_influences_marginals==1
                F = F/2;
            end

                %if the marginal can be influenced by init_bias use the Ptb and Pxt
                %calculated above, otherwise recalculate Pxt and Ptb without using the
                %init_bias
            if init_bias_influences_marginals==1,
                Pxt = Pxt_coh(:,:,1,c); % xmesh * time * stim * signed_coh
                Ptb = Ptb_coh(:,:,1,c); % time * bound * stim * signed_coh
            else
                [b_change, ~] = bcinit(B,t,dt,boundtype,uInf,tau05);
                mu = k * coh_set(c);    %no stimulation and no initial bias
                uinit = delta;          %start with delta function
                    % ?Sigma w/ coh!
                    variance = sigma^2 * (1 + coh_var_term*abs(coh_set(c)));
                [~, ~, Ptb_, ~, Pxt_] = FP4(xmesh, uinit, mu, sqrt(variance), b_change, b_margin, dt);
                    %store the arrays for each coherence level, use 1 ms time resolution 
                Ptb(:,1) = local_sum(Ptb_(2:end,1),round(1/dt));  %lower bound
                Ptb(:,2) = local_sum(Ptb_(2:end,2),round(1/dt));  %upper bound
                Pxt = Pxt_(:,1/dt:1/dt:end);
            end

                %now use Ptb and Pxt to calculate/update the marginal. 
            if coh_set(c)>0         %rightward motion
                Ptb_marginal(:,:,1) = Ptb_marginal(:,:,1) + F*Ptb;
                Pxt_marginal(:,:,1) = Pxt_marginal(:,:,1) + F*Pxt;
            elseif coh_set(c)<0     %leftward motion
                Ptb_marginal(:,:,2) = Ptb_marginal(:,:,2) + F*Ptb;
                Pxt_marginal(:,:,2) = Pxt_marginal(:,:,2) + F*Pxt;
            else                    %ambiguous motion
                Ptb_marginal(:,:,1) = Ptb_marginal(:,:,1) + 0.5*F*Ptb;
                Pxt_marginal(:,:,1) = Pxt_marginal(:,:,1) + 0.5*F*Pxt;
                Ptb_marginal(:,:,2) = Ptb_marginal(:,:,2) + 0.5*F*Ptb;
                Pxt_marginal(:,:,2) = Pxt_marginal(:,:,2) + 0.5*F*Pxt;
            end

                %repeat for the stimulated trials if needed
            if stim_bias_influences_marginals==0,
                continue;
            end;
            if init_bias_influences_marginals==1,
                Pxt = Pxt_coh(:,:,2,c);
                Ptb = Ptb_coh(:,:,2,c);
            else
                [b_change, ~] = bcinit(B,t,dt,boundtype,uInf,tau05);
                mu = k * (coh_set(c)+uStim_eff);    %no initial bias, but include stimulation 
                uinit = delta;          %start with delta function
                    % ?Sigma w/ coh!
                    variance = sigma^2 * (1 + coh_var_term*abs(coh_set(c)));
                [~, ~, Ptb_, ~, Pxt_] = FP4(xmesh, uinit, mu, sqrt(variance), b_change, b_margin, dt);
                    %store the arrays for each coherence level, use 1 ms time resolution
                Ptb(:,1) = local_sum(Ptb_(2:end,1),round(1/dt));  %lower bound
                Ptb(:,2) = local_sum(Ptb_(2:end,2),round(1/dt));  %upper bound
                Pxt = Pxt_(:,1/dt:1/dt:end);        
            end
                %now use Ptb and Pxt to calculate/update the marginal. 
            if coh_set(c)>0         %rightward motion
                Ptb_marginal(:,:,1) = Ptb_marginal(:,:,1) + F*Ptb;
                Pxt_marginal(:,:,1) = Pxt_marginal(:,:,1) + F*Pxt;
            elseif coh_set(c)<0     %leftward motion
                Ptb_marginal(:,:,2) = Ptb_marginal(:,:,2) + F*Ptb;
                Pxt_marginal(:,:,2) = Pxt_marginal(:,:,2) + F*Pxt;
            else                    %ambiguous motion
                Ptb_marginal(:,:,1) = Ptb_marginal(:,:,1) + 0.5*F*Ptb;
                Pxt_marginal(:,:,1) = Pxt_marginal(:,:,1) + 0.5*F*Pxt;
                Ptb_marginal(:,:,2) = Ptb_marginal(:,:,2) + 0.5*F*Ptb;
                Pxt_marginal(:,:,2) = Pxt_marginal(:,:,2) + 0.5*F*Pxt;
            end

        end
                
        %find out which combination of decision variable and time must be associated with
        %choosing the sure target based on logodds_crit
            %if bound crossing does not happen
        choose_strg_xt = nan(length(xmesh), max_dur);
        I = xmesh>0;
        choose_strg_xt(I,:) = log(Pxt_marginal(I,:,1)./Pxt_marginal(I,:,2)) < logodds_crit;
        I = xmesh<0;
        choose_strg_xt(I,:) = log(Pxt_marginal(I,:,2)./Pxt_marginal(I,:,1)) < logodds_crit;
        I = xmesh==0;
        if sum(I)==1;
            %choose_strg_xt(I,:) = 0<logodds_crit;
            choose_strg_xt(I,:) = mean([log(Pxt_marginal(find(I)+1,:,2)./Pxt_marginal(find(I)+1,:,1)) ; log(Pxt_marginal(find(I)-1,:,2)./Pxt_marginal(find(I)-1,:,1))]) < logodds_crit;
        end

            %if bound crossing happens 
        choose_strg_tb = nan(max_dur, 2);            
        choose_strg_tb(:,1) = abs(log(Ptb_marginal(:,1,2)./Ptb_marginal(:,1,1))) < logodds_crit;     %lower bound crossing
        choose_strg_tb(:,2) = abs(log(Ptb_marginal(:,2,1)./Ptb_marginal(:,2,2))) < logodds_crit;     %upper bound crossing

            %repeat to set up alternate criterion for stim trials (when using thetaOffset)
                %if bound crossing does not happen
        choose_strg_xt_stim = nan(length(xmesh), max_dur);
        I = xmesh>0;
        choose_strg_xt_stim(I,:) = abs(log(Pxt_marginal(I,:,1)./Pxt_marginal(I,:,2))) < logodds_crit + thetaOffset;
        I = xmesh<0;
        choose_strg_xt_stim(I,:) = abs(log(Pxt_marginal(I,:,2)./Pxt_marginal(I,:,1))) < logodds_crit + thetaOffset;
        I = xmesh==0;
        if sum(I)==1
            %choose_strg_xt_stim(I,:) = 0<logodds_crit+(thetaOffset*sign(coh_set(c)));
            choose_strg_xt_stim(I,:) = mean([abs(log(Pxt_marginal(find(I)+1,:,2)./Pxt_marginal(find(I)+1,:,1))) ; abs(log(Pxt_marginal(find(I)-1,:,2)./Pxt_marginal(find(I)-1,:,1)))]) < logodds_crit;
        end
            %if bound crossing happens 
        choose_strg_tb_stim = nan(max_dur, 2);            
        choose_strg_tb_stim(:,1) = abs(log(Ptb_marginal(:,1,2)./Ptb_marginal(:,1,1))) < logodds_crit + thetaOffset;     %lower bound crossing
        choose_strg_tb_stim(:,2) = abs(log(Ptb_marginal(:,2,1)./Ptb_marginal(:,2,2))) < logodds_crit + thetaOffset;     %upper bound crossing

        %loop through all motion strengths and trials to calculate the associated error value 
        
        minP = 1e-300;
        R = unique(trialData.CorrTrg(trialData.Coh>0));    %which choice is the rightward choice, must be 1
        % LL = 0;
        LL = zeros(size(trialData.Coh));
        expectedPright = LL;
        expectedPS = LL;
        LL2 = LL;
        t = (1 : max_dur)';
        
        try
            for s = 0 : 1,
                for c = 1 : length(coh_set),

                    Pxt = Pxt_coh(:,:,s+1,c); % xmesh * time * stim * signed_coh
                    Ptb = Ptb_coh(:,:,s+1,c); % time * bound * stim * signed_coh

                        %walk through trials one by one, due to different durations of
                        %stimulus and delays we cannot group the trials 
                    for Ts = 0 : 1,
                        trial_ind = find(trialData.uStim==s & trialData.Coh==coh_set(c) & trialData.SureT==Ts)';
                        StimDur = trialData.StimDur; Resp = trialData.Resp;

                        for i = trial_ind

                            dur = StimDur(i);
                                %calculate Pright, Pleft and Pstrg
                            if Ts == 0
                                if sum(xmesh==0)==1
                                    Pright = sum(Ptb(t<=dur,2)) + sum(Pxt(xmesh>0,dur)) + 0.5*Pxt(xmesh==0,dur);
                                    Pleft = sum(Ptb(t<=dur,1)) + sum(Pxt(xmesh<0,dur)) + 0.5*Pxt(xmesh==0,dur);
                                    Pstrg = 0;
                                else
                                    Pright = sum(Ptb(t<=dur,2)) + sum(Pxt(xmesh>0,dur));
                                    Pleft = sum(Ptb(t<=dur,1)) + sum(Pxt(xmesh<0,dur));
                                    Pstrg = 0;
                                end
                            elseif Ts == 1
                                if s==1 % allow separate Ts criteria for stim/nostim
                                    if sum(xmesh==0)==1
                                        Pright = sum(Ptb((t<=dur)&(choose_strg_tb_stim(:,2)==0),2)) + sum(Pxt((xmesh>0)&(choose_strg_xt_stim(:,dur)==0),dur)) + 0.5*sum(Pxt((xmesh==0)&(choose_strg_xt_stim(:,dur)==0),dur));
                                        Pleft = sum(Ptb((t<=dur)&(choose_strg_tb_stim(:,1)==0),1)) + sum(Pxt((xmesh<0)&(choose_strg_xt_stim(:,dur)==0),dur)) + 0.5*sum(Pxt((xmesh==0)&(choose_strg_xt_stim(:,dur)==0),dur));
                                        Pstrg = sum(Ptb((t<=dur)&(choose_strg_tb_stim(:,2)==1),2)) + sum(Ptb((t<=dur)&(choose_strg_tb_stim(:,1)==1),1)) + sum(Pxt(choose_strg_xt_stim(:,dur)==1,dur));
                                    else
                                        Pright = sum(Ptb((t<=dur)&(choose_strg_tb_stim(:,2)==0),2)) + sum(Pxt((xmesh>0)&(choose_strg_xt_stim(:,dur)==0),dur));
                                        Pleft = sum(Ptb((t<=dur)&(choose_strg_tb_stim(:,1)==0),1)) + sum(Pxt((xmesh<0)&(choose_strg_xt_stim(:,dur)==0),dur));
                                        Pstrg = sum(Ptb((t<=dur)&(choose_strg_tb_stim(:,2)==1),2)) + sum(Ptb((t<=dur)&(choose_strg_tb_stim(:,1)==1),1)) + sum(Pxt(choose_strg_xt_stim(:,dur)==1,dur));
                                    end

                                else
                                    if sum(xmesh==0)==1
                                        Pright = sum(Ptb((t<=dur)&(choose_strg_tb(:,2)==0),2)) + sum(Pxt((xmesh>0)&(choose_strg_xt(:,dur)==0),dur)) + 0.5*sum(Pxt((xmesh==0)&(choose_strg_xt(:,dur)==0),dur));
                                        Pleft = sum(Ptb((t<=dur)&(choose_strg_tb(:,1)==0),1)) + sum(Pxt((xmesh<0)&(choose_strg_xt(:,dur)==0),dur)) + 0.5*sum(Pxt((xmesh==0)&(choose_strg_xt(:,dur)==0),dur));
                                        Pstrg = sum(Ptb((t<=dur)&(choose_strg_tb(:,2)==1),2)) + sum(Ptb((t<=dur)&(choose_strg_tb(:,1)==1),1)) + sum(Pxt(choose_strg_xt(:,dur)==1,dur));
                                    else
                                        Pright = sum(Ptb((t<=dur)&(choose_strg_tb(:,2)==0),2)) + sum(Pxt((xmesh>0)&(choose_strg_xt(:,dur)==0),dur));
                                        Pleft = sum(Ptb((t<=dur)&(choose_strg_tb(:,1)==0),1)) + sum(Pxt((xmesh<0)&(choose_strg_xt(:,dur)==0),dur));
                                        Pstrg = sum(Ptb((t<=dur)&(choose_strg_tb(:,2)==1),2)) + sum(Ptb((t<=dur)&(choose_strg_tb(:,1)==1),1)) + sum(Pxt(choose_strg_xt(:,dur)==1,dur));
                                    end
                                end
                            end

                            Ptot = Pright + Pleft + Pstrg;
                            if abs(Ptot-1)>1e-3
                                warning('FitDiffusion:InadequatePropagation', 'the probabilities (Pright, Pleft, and Pstrg) do not add up to one!\nlonger propagation is needed\n\tPtot=%f (for c: %d, Ts: %d, n: %d)\n', Ptot, trialData.Coh(n), trialData.SureT(n), n);
                            end
                            Pright = Pright/Ptot;
                            Pleft = Pleft/Ptot;
                            Pstrg = Pstrg/Ptot;

                            %calculate the expected probability correct and probability sure target 
                            expectedPright(i) = Pright/(Pright+Pleft);
                            expectedPS(i) = Pstrg;

                            % "secondary error" -- not what is being minimized, but something to
                            % measure for questions of 'goodness of prediction' (model comparison)
                            if s==1 % for now record err of Psure on stim trials only
                                LL2(i) = (Resp(i)==3)*log(max(Pstrg,minP)) + (Resp(i)==R||Resp(i)==3-R)*log(max(Pright+Pleft,minP));
                            end
                            
                            switch fitwhat{s+1}
                                case 'PrightPS'
                                        % fit Pright and PS
                                    LL(i) = (Resp(i)==3)*log(max(Pstrg,minP)) + (Resp(i)==R)*log(max(Pright,minP)) + (Resp(i)==3-R)*log(max(Pleft,minP));
                                case 'Pright'
                                        % fit only Pright. Renormalize Pright for trials in which Ts was shown
                                    Pright = Pright/(Pright+Pleft);
                                    LL(i) = (Resp(i)==R)*log(max(Pright,minP)) + (Resp(i)==3-R)*log(max(1-Pright,minP));
                                otherwise
                                    error('unknown fitwhat!');
                            end

                        end

                    end
                end
            end
        catch E
            fprintf('\n\nERROR!!!\n');
            fprintf('call %d, param %g %g %g %g %g %g\n', callNum_, k, B, logodds_crit, weber_fraction, uStim_eff, init_bias);
        %     fprintf('c: %d, Ts: %d, uStim: %d, i: %d\n', c, Ts, s, i);
            rethrow(E);
        end
        
        err = -sum(LL);
        err2 = -sum(LL2);
        n = length(LL);
        
        % print progress report!
        fprintf ( '\n\n\n****************************************\n');
        fprintf ( 'file: %s\n' , filename );
        fprintf ( 'run %d\n' , callNum_ );
        fprintf ( '\tk= %f,\n\tsigma= %f\n\tlogodds_crit= %f,\n\tuStim_effect= %f\n\tinit_bias= %f\n\tthetaOffset= %f\n\tkmult= %f\n\tBmult= %f\n\tvarOffset= %f\n\tCohVar= %f\n\tCVoffset= %f\n', k, sigma, logodds_crit, uStim_eff, init_bias, thetaOffset, stim_k_mult, stim_B_mult, varOffset, coh_var_term, coh_var_offset);
        fprintf ( 'err: %f\n' , err );
        fprintf ('number of processed trials: %d\n', n);
        
        callNum_ = callNum_ + 1;
    end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function param2 = getParam(param1,~,~)
        if ~all(fixed==1)
            param2(fixed==0) = param1;              %get adjustable parameters from param1
        end
        param2(fixed==1) = guess(fixed==1);     %get fixed parameters from guess
    end


end
