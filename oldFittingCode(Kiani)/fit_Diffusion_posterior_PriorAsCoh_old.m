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
        
        tic
        
        spect = 0;
        boundtype = 'flat';
%         boundtype = 'hyperbolic';       
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The model has the following free parameters
        %   K   the factor that converts the motion strength to evidence   e=kC
        %   B   the bound height
        %   logodds_crit
        %   weber_fraction
        %   uStim_eff
        %   init_bias
        %   thetaOffset
        %   stim_choice_bias (scb)
        %   k_mult
        %   B_mult
        
        init_bias_influences_marginals = 0;
        stim_bias_influences_marginals = 0;

        param = getParam(param);

        % % param scaling (see calling function)
        % param(2) = param(2)*100;
        % param(3) = param(3)*2;
        % param(6) = param(6)/10;
                    % abandoned: did not seem to matter -CF

        k = param(1);
        B = abs(param(2));  %don't accept negative bound heights
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
        if spect && varOffset~=1
            error('spectral method cannot currently handle varOffset!');
        end
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
        if numel(param)>13
            x0 = param(14); % starting point (as fraction of bound)
        else
            x0 = 0;
        end
        % done w/ setting free params
        %***************************************************
        %***************************************************
        
        sigma = 1;
        
                %define the grid resolution for time and decision variable
        if spect
            dt = 0.5; % in milliseconds!
        else
            dt = 0.1;
        end
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
%         delta(abs(xmesh)==min(abs(xmesh))) = 1;
            %NEW: adjust for x0
        delta(abs(xmesh-x0*B) == min(abs(xmesh-x0*B))) = 1;

            %find the probability density of decision variable across time for stimulated and
            %nonstimulated trials 
        Ptb_coh = zeros(max_dur, 2, 2, length(coh_set));                % time * bound * stim * signed_coh
        Pxt_coh = zeros(length(xmesh), max_dur, 2, length(coh_set));    % xmesh * time * stim * signed_coh

        % make versions of xmesh etc for spectral_dtb (length must be power of 2, and B is off by factor of sqrt(1000))
        if spect
            kS = k*sqrt(1000); % convert for spect/analytic, which uses seconds instead of ms
            B_S = B/sqrt(1000);
            % xmeshS = ((xmesh(1):(xmesh(end)-xmesh(1))/511:xmesh(end))/100)';
            drift = [kS*(coh_set+0*uStim_eff) + init_bias ; kS*(coh_set+1*uStim_eff) + init_bias];
            dtS = dt/1000; % change to seconds
            tS = t/1000;
            xmesh = linspace(-B_S-max(abs(drift))*dtS-4*sqrt(dtS), B_S+max(abs(drift))*dtS+4*sqrt(dtS), 512)';
            deltaS = zeros(size(xmesh));
%             deltaS(abs(xmesh)==min(abs(xmesh))) = 1;
                %NEW: adjust for x0
            deltaS(abs(xmesh-x0*B) == min(abs(xmesh-x0*B))) = 1;
            Ptb_S_coh = zeros(max_dur, 2, 2, length(coh_set));                % time * bound * stim * signed_coh
            Pxt_S_coh = zeros(length(xmesh), max_dur, 2, length(coh_set));    % xmesh * time * stim * signed_coh
        end

        % FP4elapsed = 0; spectElapsed = 0; anaElapsed = 0;
        for s = 0 : 1
                %do not run FP4 unless it is really necessary 
            if stim_bias_influences_marginals==0
                if init_bias_influences_marginals==0 && s==0 && all(trialData.uStim==1)
                    continue;
                elseif s==1 && all(trialData.uStim==0)
                    continue;
                end
            end
            if ~spect
                %run FP4 to calculate Ptb and Pxt for each coherence and stimulation condition 
            for c = 1 : length(coh_set)
                if s==1 && uStim_eff==0     %do not need to recalculate if uStim_eff is zero
                    Ptb_coh(:,:,s+1,c) = Ptb_coh(:,:,1,c);
                    Pxt_coh(:,:,s+1,c) = Pxt_coh(:,:,1,c);
                else
                    % must re-initialize b_change, which (ahem) changes after each call of FP4
                    [b_change, ~] = bcinit(B,t,dt,boundtype,uInf,tau05);
                    if s==1
                        mu = k * stim_k_mult * (coh_set(c)+s*uStim_eff) + init_bias;
                        b_change(:,1) = b_change(:,1)*min(stim_B_mult,1) - (min(stim_B_mult,1)*B-B); % implements stim_B_mult using
                        b_change(:,2) = b_change(:,2)*min(stim_B_mult,1) + (min(stim_B_mult,1)*B-B); % existing b_change functionality
                    else
                        mu = k * (coh_set(c)+s*uStim_eff) + init_bias;
                    end
                    uinit = delta;
    %                 FP4start = tic;
    
%                     [~, ~, Ptb_, ~, Pxt_] = FP4(xmesh, uinit, mu, sigma+s*(varOffset*sigma-sigma), b_change, b_margin, dt);
%                     [~, ~, Ptb_, ~, Pxt_] = FP4(xmesh, uinit, mu, sqrt(sigma^2+s*varOffset), b_change, b_margin, dt);

%                         coh-dependent sigma!
%                       variance = (sigma^2+varOffset*s) * (1 + (coh_var_term+coh_var_offset*s)*abs(coh_set(c)))
                      variance = (sigma^2+varOffset*s) + coh_var_term*abs(coh_set(c));

                      [~, ~, Ptb_, ~, Pxt_] = FP4(xmesh, uinit, mu, sqrt(variance), b_change, b_margin, dt);
%                 FP4elapsed = FP4elapsed + toc(FP4start);
                        %store the arrays for each coherence level, use 1 ms time resolution
                    Ptb_coh(:,1,s+1,c) = local_sum(Ptb_(2:end,1),1/dt);  %lower bound
                    Ptb_coh(:,2,s+1,c) = local_sum(Ptb_(2:end,2),1/dt);  %upper bound
                    Pxt_coh(:,:,s+1,c) = Pxt_(:,1/dt:1/dt:end);
                end
            end

            else
%             spectStart = tic;
            [~, B_S] = bcinit(B,t,dt,boundtype,uInf,tau05);
            if s==1
                mu = kS * stim_k_mult * (coh_set+s*uStim_eff) + init_bias;
                B_S = B_S * min(stim_B_mult,1);
            else
                mu = kS * (coh_set+s*uStim_eff) + init_bias;
            end
            D = spectral_dtbInd(mu', (dtS:dtS:(length(tS)*dtS))', B_S, -B_S, xmesh, deltaS, 1);
            Ptb_S = [D.lo.pdf_t D.up.pdf_t];
            Pxt_S = D.notabs.pdf;
            Ptb_S_coh(:,1,s+1,:) = local_sum(Ptb_S(:,1:length(coh_set)),0.001/dtS);  %lower bound
            Ptb_S_coh(:,2,s+1,:) = local_sum(Ptb_S(:,length(coh_set)+1:end),0.001/dtS);  %upper bound
            Pxt_S_coh(:,:,s+1,:) = permute(Pxt_S(:,:,0.001/dtS:0.001/dtS:end),[2 3 1]);
%             spectElapsed = spectElapsed + toc(spectStart)
            end

        end

        if spect
            Ptb_coh = Ptb_S_coh;
            Pxt_coh = Pxt_S_coh;
        end
        
        
%         cohind = 22;
%         % take a look at prob density across time for coh=coh_set(cohind) and:
%         figure; set(gcf,'Position',[200 600 500 800]);
%         temp = squeeze(Pxt_coh(:,:,1,cohind)); % nostim
%         temp(temp<1e-10) = 1e-10;
%         subplot(3,1,1); plot(Ptb_coh(:,2,1,cohind)); ylim([0 1.5e-3]);
%         subplot(3,1,2); contourf(log10(temp)); title('log10(PxtFull), coh=0, nostim'); colorbar;
%         subplot(3,1,3); plot(Ptb_coh(:,1,1,cohind)); ylim([0 1.5e-3]);
%         figure; set(gcf,'Position',[500 600 500 800]);
%         temp = squeeze(Pxt_coh(:,:,2,cohind)); % stim
%         temp(temp<1e-10) = 1e-10;
%         subplot(3,1,1); plot(Ptb_coh(:,2,2,cohind)); ylim([0 1.5e-3]);
%         subplot(3,1,2); contourf(log10(temp)); title('log10(PxtFull), coh=0, stim'); colorbar;
%         subplot(3,1,3); plot(Ptb_coh(:,1,2,cohind)); ylim([0 1.5e-3]);
%         
%         pause;

%         % sanity check for marginals
%         for c = 1:2:length(coh_set)
%             temp = squeeze(Pxt_coh(:,:,1,c)); % nostim
%             temp(temp<1e-10) = 1e-10;
%             figure(5423);clf; contourf(log10(temp)+10); title(sprintf('PxtFull,coh=%d,nostim',coh_set(c))); colorbar;
%             pause;
%         end
%         for c = 1:2:length(coh_set)
%             temp = squeeze(Pxt_coh(:,:,2,c)); % stim
%             temp(temp<1e-10) = 1e-10;
%             figure(5423);clf; contourf(log10(temp)+10); title(sprintf('PxtFull,coh=%d,stim',coh_set(c))); colorbar;
%             pause;
%         end

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
                [b_change, B_S] = bcinit(B,t,dt,boundtype,uInf,tau05);
                if ~spect
                    mu = k * coh_set(c);    %no stimulation and no initial bias
                    delta = zeros(size(xmesh));
                    delta(abs(xmesh)==min(abs(xmesh))) = 1;
                    uinit = delta;          %start with delta function (no x0)
                        % ?Sigma w/ coh!
                        variance = sigma^2 * (1 + coh_var_term*abs(coh_set(c)));
                    [~, ~, Ptb_, ~, Pxt_] = FP4(xmesh, uinit, mu, sqrt(variance), b_change, b_margin, dt);
                        %store the arrays for each coherence level, use 1 ms time resolution 
                    Ptb(:,1) = local_sum(Ptb_(2:end,1),round(1/dt));  %lower bound
                    Ptb(:,2) = local_sum(Ptb_(2:end,2),round(1/dt));  %upper bound
                    Pxt = Pxt_(:,1/dt:1/dt:end);
                else
                    mu = kS * coh_set(c);
                    deltaS = zeros(size(xmesh));
                    deltaS(abs(xmesh)==min(abs(xmesh))) = 1; %start with delta function (no x0)
                    D = spectral_dtbInd(mu', (dtS:dtS:(length(t)*dtS))', B_S, -B_S, xmesh, deltaS, 1);
                    Ptb_ = [D.lo.pdf_t D.up.pdf_t];
                    Pxt_ = squeeze(D.notabs.pdf);
                        %store the arrays for each coherence level, use 1 ms time resolution 
                    Ptb(:,1) = local_sum(Ptb_(:,1),0.001/dtS);  %lower bound
                    Ptb(:,2) = local_sum(Ptb_(:,2),0.001/dtS);  %upper bound
                    Pxt = Pxt_(:,0.001/dtS:0.001/dtS:end);
                end
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
                [b_change, B_S] = bcinit(B,t,dt,boundtype,uInf,tau05);
                if ~spect
                    mu = k * (coh_set(c)+uStim_eff);    %no initial bias, but include stimulation 
                    uinit = delta;          %start with delta function (no x0)
                        % ?Sigma w/ coh!
                        variance = sigma^2 * (1 + coh_var_term*abs(coh_set(c)));
                    [~, ~, Ptb_, ~, Pxt_] = FP4(xmesh, uinit, mu, sqrt(variance), b_change, b_margin, dt);
                else
                    mu = kS * (coh_set(c)+uStim_eff);    %no initial bias, but include stimulation 
                    D = spectral_dtbInd(mu', (dtS:dtS:(length(t)*dtS))', B_S, -B_S, xmesh, deltaS, 1);
                    Ptb_ = [D.lo.pdf_t D.up.pdf_t];
                    Pxt_ = squeeze(D.notabs.pdf);
                end                
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

%         temp = Pxt_marginal(:,:,1);
%         temp(temp<1e-10) = 1e-10;
%         figure(123); subplot(3,2,1); contourf(log10(temp)); title('log(PxtMarg), rightward'); colorbar;
%         temp = Pxt_marginal(:,:,2);
%         temp(temp<1e-10) = 1e-10;
%         figure(123); subplot(3,2,2); contourf(log10(temp)); title('log(PxtMarg), leftward'); colorbar;

        
%         % force left and right marginals to have equal probability for a given time slice -CRF
          % codebank(39)
          
%         %smooth the probability densitiy of decision variable based on temporal
%         %uncertainty.  -- codebank(40)

        % for reference:
        %         Ptb_marginal: time * bound * motion_direction (1 right, 2 left)
        %         Pxt_marginal: xmesh * time * motion_direction
                
                    %find out which combination of decision variable and time must be associated with
                    %choosing the sure target based on logodds_crit
                        %if bound crossing does not happen
                choose_strg_xt = nan(length(xmesh), max_dur);
                I = xmesh>0;
                choose_strg_xt(I,:) = log(Pxt_marginal(I,:,1)./Pxt_marginal(I,:,2)) < logodds_crit;
%                     asdf1 = log(Pxt_marginal(I,:,1)./Pxt_marginal(I,:,2));        
                I = xmesh<0;
                choose_strg_xt(I,:) = log(Pxt_marginal(I,:,2)./Pxt_marginal(I,:,1)) < logodds_crit;
%                     asdf2 = log(Pxt_marginal(I,:,2)./Pxt_marginal(I,:,1)); 
                I = xmesh==0;
                if sum(I)==1;
                    %choose_strg_xt(I,:) = 0<logodds_crit;
                    choose_strg_xt(I,:) = mean([log(Pxt_marginal(find(I)+1,:,2)./Pxt_marginal(find(I)+1,:,1)) ; log(Pxt_marginal(find(I)-1,:,2)./Pxt_marginal(find(I)-1,:,1))]) < logodds_crit;
%                         asdf3 = mean([log(Pxt_marginal(find(I)+1,:,2)./Pxt_marginal(find(I)+1,:,1)) ; log(Pxt_marginal(find(I)-1,:,2)./Pxt_marginal(find(I)-1,:,1))]); 
%                     asdf = [asdf2;asdf3;asdf1];
%                 else
%                     asdf = [asdf2;asdf1];
                end
%                     figure(123); subplot(3,2,3); contourf(repmat(1:size(asdf,2),size(asdf,1),1)',repmat(xmesh,1,size(asdf,2))',asdf'); colorbar; title('Log odds correct');
%                     figure(123); subplot(3,2,4); contourf(repmat(1:size(choose_strg_xt,2),size(choose_strg_xt,1),1)',repmat(xmesh,1,size(choose_strg_xt,2))',choose_strg_xt'); colorbar; title('choose Ts or not');

                    %if bound crossing happens 
                choose_strg_tb = nan(max_dur, 2);            
                choose_strg_tb(:,1) = abs(log(Ptb_marginal(:,1,2)./Ptb_marginal(:,1,1))) < logodds_crit;     %lower bound crossing
                choose_strg_tb(:,2) = abs(log(Ptb_marginal(:,2,1)./Ptb_marginal(:,2,2))) < logodds_crit;     %upper bound crossing

                    %repeat to set up alternate criterion for stim trials (when using thetaOffset)
                        %if bound crossing does not happen
                choose_strg_xt_stim = nan(length(xmesh), max_dur);
                I = xmesh>0;
                choose_strg_xt_stim(I,:) = abs(log(Pxt_marginal(I,:,1)./Pxt_marginal(I,:,2))) < logodds_crit + thetaOffset;
%                     asdf1 = log(Pxt_marginal(I,:,1)./Pxt_marginal(I,:,2));        
                I = xmesh<0;
                choose_strg_xt_stim(I,:) = abs(log(Pxt_marginal(I,:,2)./Pxt_marginal(I,:,1))) < logodds_crit - thetaOffset;
%                     asdf2 = log(Pxt_marginal(I,:,2)./Pxt_marginal(I,:,1)); 
                I = xmesh==0;
                if sum(I)==1
                    %choose_strg_xt_stim(I,:) = 0<logodds_crit+(thetaOffset*sign(coh_set(c)));
                    choose_strg_xt_stim(I,:) = mean([abs(log(Pxt_marginal(find(I)+1,:,2)./Pxt_marginal(find(I)+1,:,1))) ; abs(log(Pxt_marginal(find(I)-1,:,2)./Pxt_marginal(find(I)-1,:,1)))]) < logodds_crit;
%                         asdf3 = mean([log(Pxt_marginal(find(I)+1,:,2)./Pxt_marginal(find(I)+1,:,1)) ; log(Pxt_marginal(find(I)-1,:,2)./Pxt_marginal(find(I)-1,:,1))]); 
%                     asdf = [asdf2;asdf3;asdf1];
%                 else
%                     asdf = [asdf2;asdf1];
                end
%                     figure(123); subplot(3,2,5); contourf(repmat(1:size(asdf,2),size(asdf,1),1)',repmat(xmesh,1,size(asdf,2))',asdf'); colorbar; title('Log odds correct STIM');
%                     figure(123); subplot(3,2,6); contourf(repmat(1:size(choose_strg_xt_stim,2),size(choose_strg_xt_stim,1),1)',repmat(xmesh,1,size(choose_strg_xt_stim,2))',choose_strg_xt_stim'); colorbar; title('choose Ts or not STIM');

                    %if bound crossing happens 
                choose_strg_tb_stim = nan(max_dur, 2);            
                choose_strg_tb_stim(:,1) = abs(log(Ptb_marginal(:,1,2)./Ptb_marginal(:,1,1))) < logodds_crit + thetaOffset;     %lower bound crossing
                choose_strg_tb_stim(:,2) = abs(log(Ptb_marginal(:,2,1)./Ptb_marginal(:,2,2))) < logodds_crit - thetaOffset;     %upper bound crossing

        %loop through all motion strengths and trials to calculate the associated error value 
        
        minP = 1e-300;
        R = unique(trialData.CorrTrg(trialData.Coh>0));    %which choice is the rightward choice, must be 1
        % LL = 0;
        LL = zeros(size(trialData.Coh));
        expectedPright = LL;
        expectedPS = LL;
        LL2 = LL;
        t = (1 : max_dur)';
        
%         tic
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
%         toc
        
        err = -sum(LL);
        err2 = -sum(LL2);
        n = length(LL);
        
        % print progress report!
        fprintf ( '\n\n\n****************************************\n');
        fprintf ( 'file: %s\n' , filename );
        fprintf ( 'run %d\n' , callNum_ );
        fprintf ( '\tk= %f,\n\tB= %f\n\tlogodds_crit= %f,\n\tuStim_effect= %f\n\tinit_bias= %f\n\tthetaOffset= %f\n\tkmult= %f\n\tBmult= %f\n\tvarOffset= %f\n\tCohVar= %f\n\tCVoffset= %f\n\tx0= %f\n', k, B, logodds_crit, uStim_eff, init_bias, thetaOffset, stim_k_mult, stim_B_mult, varOffset, coh_var_term, coh_var_offset, x0);
        fprintf ( 'err: %f\n' , err );
        fprintf ('number of processed trials: %d\n', n);
        
        callNum_ = callNum_ + 1;
        toc
    end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function param2 = getParam(param1,~,~)
        if ~all(fixed==1)
            param2(fixed==0) = param1;              %get adjustable parameters from param1
        end
        param2(fixed==1) = guess(fixed==1);     %get fixed parameters from guess
    end


end
