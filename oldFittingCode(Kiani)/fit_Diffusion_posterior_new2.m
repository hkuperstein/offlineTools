function [model_param, modelLL, trial_data] = fit_Diffusion_posterior_new2(trial_data, guess, fixed, feedback, options)

%     default.fitwhat = {'PrightPS', 'PrightPS'};
    default.fitwhat = {'PrightPS', 'Pright'};
    default.init_bias_influences_marginals = false;
    default.stim_bias_influences_marginals = false;
    default.use_midline_as_criterion = true;   % set to 1 to use conventional mapping: positive decision variable means right choice and negative means left choice
                                                % set to 0 to depart from the conventional mapping and instead choose based on the sign of logodds_right 
    default.parfor = true;
    default.plot = false;
    default.coh_set = [];
    default.coh_set_freq = [];
    if nargin<5 || isempty(options)
        opt = default;
    else
        opt = safeStructAssign(default, options);
    end
    
    if nargin<4 || isempty(feedback)
        feedback = 1;
    end
    
    call_num = 1;
    model_param = struct('init', guess, 'fixed', fixed, 'final', [], 'se', []);
    if all(fixed==1)
        fh = -1;
        fprintf('\n\nall parameters are fixed, no optimization will be performed\n\n');
        modelLL = -fitmodel_MLEerr([]);
        return;
    end
    
    if feedback>0,
        fh = figure;
        hold on;   
        xlabel('Call number' );
        ylabel('-LL');
    end
    
    options = optimset('Display' ,'final' ,'MaxFunEvals' ,500*sum(fixed==0) ,'MaxIter' ,500*sum(fixed==0), 'TolX', 1e-3, 'TolFun', 1e-2);
    [p, fval] = fminsearch(@fitmodel_MLEerr, guess(fixed==0), options);
    model_param.final = getParam(p, guess, fixed);
    modelLL = -fval;
    
    feedback = 0;
    Hessian = calcHessian(@fitmodel_MLEerr, model_param.final(model_param.fixed~=1));
    model_param.hessian = Hessian;
    model_param.se = nan(size(model_param.final));
    model_param.se(model_param.fixed~=1) = diag(sqrt(inv(Hessian)));
    
    guess = model_param.final;
    fixed = ones(size(guess));
    opt.plot = 0;
    fitmodel_MLEerr([]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The model has the following free parameters
    %   K   the factor that converts the motion strength to evidence   e=kC
    %   B   the bound height
    %   theta
    %   waber_fraction
    %   uStim_eff
    %   
        %this diagram shows the events in a trials based on the model and tells how reactions times should be calculated     
        %
        %           Dots On                                    Dots Off     Go  
        %               |           StimDur                        |  Delay |  T1 (detection of Go signal) 
        %               |------------------------------------------|........|........|
        %
        %                               Diffusion Start                             
        %                   T0                |     StimDur=Maximum diffusion duration
        %               |.....................|------------------------------------------|
        %   
        %                                                                                T2 (Motor Delay) 
        %                                                                              |...........|
        % Reaction time is calculated from Dots Off (go signal)
        % T2 is modified based on the bound crossing time
        %
    
    function err = fitmodel_MLEerr(param)
        param = getParam(param,guess,fixed);
        
        k = param(1);
        B = abs(param(2));
        sigma = 1;
        theta = abs(param(3));
        dcoh_stim = param(4);
        dcoh_bias = param(5);
        dk = param(6);
        dsigma = param(7);
        sigma2_beta = param(8);
%         theta2 = abs(param(9));
            % for now, manually yoke theta and theta2
        theta2 = theta;
        dtheta = param(10);
        dtheta2 = param(11);
        weber_fraction = 0;
        
            %get the coherence levels, the frequency of each coherence level and the maximum
            %stimulus duration 
            % Note that coh must be signed for the algorithm to work properly 
        coh_set = unique(trial_data.Coh);
        coh_set_freq = nan(length(coh_set), 2);
        for s = 0 : 1
            for c = 1 : length(coh_set)
                coh_set_freq(c,s+1) = sum(trial_data.Coh==coh_set(c)&trial_data.uStim==s)/length(trial_data.Coh);
            end
        end
        
            %define the grid resolution for time and decision variable 
        dt = 0.1;
        dx = min(0.1, B/100);
            %define the time axis 
        max_dur = max(trial_data.StimDur);
        t = dt:dt:max_dur;
            %define xmesh and the ghost cells for bound crossing probability
        b_margin = repmat(4*max(sigma,sigma+dsigma), [1 2]);
        xmesh = (-B-b_margin(1)+dx : dx : B+b_margin(2)-dx)';
            % adjust dx so that the mesh has zero and is therefore symmetric 
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
        Ptb_coh_formarg = zeros(max_dur, 2, 2, length(coh_set));                % time * bound * stim * signed_coh
        Pxt_coh_formarg = zeros(length(xmesh), max_dur, 2, length(coh_set));    % xmesh * time * stim * signed_coh
        Ptb_marginal = zeros(max_dur, 2, 2);                            % time * bound * motion_direction(1 right, 2 left)
        Pxt_marginal = zeros(length(xmesh), max_dur, 2);                % xmesh * time * motion_direction
        for s = 0 : 1
                %run FP4 to calculate Ptb and Pxt for each coherence and stimulation condition 
            if opt.parfor
                s_parfor = s + 1;
                parfor ci = 1 : length(coh_set)
                    if coh_set_freq(ci,s_parfor)==0
                        continue;
                    end
                    mu = (k+s*dk) * (coh_set(ci)+s*dcoh_stim) + dcoh_bias;
                    uinit = delta;      %start with delta function
                    b_change = zeros(length(t),2);
                    [~, ~, Ptb, ~, Pxt] = FP4(xmesh, uinit, mu, sigma+s*dsigma+sigma2_beta*abs(coh_set(ci)), b_change, b_margin, dt);
                        %store the arrays for each coherence level, use 1 ms time resolution 
                    Ptb_coh(:,:,s_parfor,ci) = [local_sum(Ptb(2:end,1),round(1/dt)), local_sum(Ptb(2:end,2),round(1/dt))];
                    Pxt_coh(:,:,s_parfor,ci) = Pxt(:,1/dt:1/dt:end);
                end
                parfor ci = 1 : length(coh_set)
                    if coh_set_freq(ci,s_parfor)==0
                        continue;
                    end
                        %calculate the marginals
                    if (s==0 && opt.init_bias_influences_marginals==false && dcoh_bias~=0) || ...
                       (s==1 && opt.stim_bias_influences_marginals==false && opt.init_bias_influences_marginals==false)
                        mu = k * coh_set(ci);
                        uinit = delta;      %start with delta function
                        b_change = zeros(length(t),2);
                        [~, ~, Ptb, ~, Pxt] = FP4(xmesh, uinit, mu, sigma+sigma2_beta*abs(coh_set(ci)), b_change, b_margin, dt);
                        Ptb_coh_formarg(:,:,s_parfor,ci) = [local_sum(Ptb(2:end,1),round(1/dt)), local_sum(Ptb(2:end,2),round(1/dt))];
                        Pxt_coh_formarg(:,:,s_parfor,ci) = Pxt(:,1/dt:1/dt:end);
                    elseif s==1 && opt.stim_bias_influences_marginals==false && opt.init_bias_influences_marginals==true
                        Ptb_coh_formarg(:,:,s_parfor,ci) = Ptb_coh(:,:,s,ci);
                        Pxt_coh_formarg(:,:,s_parfor,ci) = Pxt_coh(:,:,s,ci);
                    else
                        Ptb_coh_formarg(:,:,s_parfor,ci) = Ptb_coh(:,:,s_parfor,ci); %#ok<*PFBNS>
                        Pxt_coh_formarg(:,:,s_parfor,ci) = Pxt_coh(:,:,s_parfor,ci);
                    end
                end
            else
                    %calculate Pxt and Ptb for the stimulated and non-stimulated trials
                for c = 1 : length(coh_set)
                    if coh_set_freq(c,s+1)==0
                        continue;
                    end
                    mu = (k+s*dk) * (coh_set(c)+s*dcoh_stim) + dcoh_bias;
                    uinit = delta;      %start with delta function
                    b_change = zeros(length(t),2);
                    [~, ~, Ptb, ~, Pxt] = FP4(xmesh, uinit, mu, sigma+s*dsigma+sigma2_beta*abs(coh_set(ci)), b_change, b_margin, dt);
                        %store the arrays for each coherence level, use 1 ms time resolution 
                    Ptb_coh(:,:,s+1,c) = [local_sum(Ptb(2:end,1),round(1/dt)), local_sum(Ptb(2:end,2),round(1/dt))];
                    Pxt_coh(:,:,s+1,c) = Pxt(:,1/dt:1/dt:end);
                end
                    %calculate Pxt and Ptb for the calculation of marginal distributions
                for c = 1 : length(coh_set)
                    if coh_set_freq(c,s+1)==0
                        continue;
                    end
                    if (s==0 && opt.init_bias_influences_marginals==false && dcoh_bias~=0) || ...
                       (s==1 && opt.stim_bias_influences_marginals==false && opt.init_bias_influences_marginals==false)
                        mu = k * coh_set(c);
                        uinit = delta;      %start with delta function
                        b_change = zeros(length(t),2);
                        [~, ~, Ptb, ~, Pxt] = FP4(xmesh, uinit, mu, sigma+sigma2_beta*abs(coh_set(ci)), b_change, b_margin, dt);
                        Ptb_coh_formarg(:,:,s+1,c) = [local_sum(Ptb(2:end,1),round(1/dt)), local_sum(Ptb(2:end,2),round(1/dt))];
                        Pxt_coh_formarg(:,:,s+1,c) = Pxt(:,1/dt:1/dt:end);
                    elseif s==1 && opt.stim_bias_influences_marginals==false && opt.init_bias_influences_marginals==true
                        Ptb_coh_formarg(:,:,s+1,c) = Ptb_coh(:,:,s,c);
                        Pxt_coh_formarg(:,:,s+1,c) = Pxt_coh(:,:,s,c);
                    else
                        Ptb_coh_formarg(:,:,s+1,c) = Ptb_coh(:,:,s+1,c);
                        Pxt_coh_formarg(:,:,s+1,c) = Pxt_coh(:,:,s+1,c);
                    end
                end
            end
                %calculate the marginals
            if s==1 && opt.stim_bias_influences_marginals==false
                continue;
            else
                for c = 1 : length(coh_set)
                    if isempty(opt.coh_set)
                        if opt.stim_bias_influences_marginals==true
                            F = coh_set_freq(c,s+1)/sum(coh_set_freq(:));
                        else
                            F = coh_set_freq(c,s+1)/sum(coh_set_freq(:,s+1));
                        end
                    else
                        if all(coh_set(c)~=opt.coh_set)
                            continue;
                        else
                            if opt.stim_bias_influences_marginals==true
                                F = opt.coh_set_freq(opt.coh_set==coh_set(c),s+1)/sum(opt.coh_set_freq(:));
                            else
                                F = opt.coh_set_freq(opt.coh_set==coh_set(c),s+1)/sum(opt.coh_set_freq(:,s+1));
                            end
                        end
                    end
                        %now use Ptb and Pxt to calculate/update the marginal. 
                    Pxt = Pxt_coh_formarg(:,:,s+1,c);
                    Ptb = Ptb_coh_formarg(:,:,s+1,c);
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
            end
        end
        
                %smooth the probability densitiy of decision variable based on temporal uncertainty. 
                %ideally the application of weber fraction should happen before reduction of the
                %temporal resolution of Pxt_marginal and Ptb_marginal, but this process is extremely
                %memory expensive, so I am settling for the 1 ms approximation.
        if weber_fraction > 0
            S = Ptb_marginal(:,:,1);
            Ptb_marginal(:,:,1) = applyWeberTime(S, weber_fraction, 0, 1);
            S = Ptb_marginal(:,:,2);
            Ptb_marginal(:,:,2) = applyWeberTime(S, weber_fraction, 0, 1);
            S = Pxt_marginal(:,:,1);
            Pxt_marginal(:,:,1) = applyWeberTime(S', weber_fraction, 0, 0)';
            S = Pxt_marginal(:,:,2);
            Pxt_marginal(:,:,2) = applyWeberTime(S', weber_fraction, 0, 0)';
        end
        
        
        
                %find out which combination of decision variable and time must be associated with
                %choosing the sure target based on theta
        if opt.use_midline_as_criterion==false
                    %if bound crossing does not happen
            LogPostOdds_xt = log(Pxt_marginal(:,:,1)./Pxt_marginal(:,:,2));
            choose_strg_xt = abs(LogPostOdds_xt) < theta;
                    %if bound crossing happens 
            LogPostOdds_tb = [log(Ptb_marginal(:,1,1)./Ptb_marginal(:,1,2)), log(Ptb_marginal(:,2,1)./Ptb_marginal(:,2,2))];    %[lower_bound, upper_bound]
            choose_strg_tb = abs(LogPostOdds_tb) < theta;
        else
                    %if bound crossing does not happen
            choose_strg_xt = nan(length(xmesh), max_dur);
            I = xmesh>0;
            choose_strg_xt(I,:) = log(Pxt_marginal(I,:,1)./Pxt_marginal(I,:,2)) < theta;
            I = xmesh<0;
            choose_strg_xt(I,:) = log(Pxt_marginal(I,:,2)./Pxt_marginal(I,:,1)) < theta2;
            I = xmesh==0;
            choose_strg_xt(I,:) = 0<theta;
                    %if bound crossing happens 
            choose_strg_tb = nan(max_dur, 2);            
            choose_strg_tb(:,1) = log(Ptb_marginal(:,1,2)./Ptb_marginal(:,1,1)) < theta;     %lower bound crossing
            choose_strg_tb(:,2) = log(Ptb_marginal(:,2,1)./Ptb_marginal(:,2,2)) < theta2;     %upper bound crossing
        end
        

            %*** MY VERSION ***%
                    %find out which combination of decision variable and time must be associated with
                    %choosing the sure target based on theta
                        %if bound crossing does not happen
                choose_strg_xt = nan(length(xmesh), max_dur);
                I = xmesh>0;
                choose_strg_xt(I,:) = log(Pxt_marginal(I,:,1)./Pxt_marginal(I,:,2)) < theta;
                I = xmesh<0;
                choose_strg_xt(I,:) = log(Pxt_marginal(I,:,2)./Pxt_marginal(I,:,1)) < theta2;
                I = xmesh==0;
                if sum(I)==1;
                    choose_strg_xt(I,:) = mean([log(Pxt_marginal(find(I)+1,:,2)./Pxt_marginal(find(I)+1,:,1)) ; log(Pxt_marginal(find(I)-1,:,2)./Pxt_marginal(find(I)-1,:,1))]) < theta;
                end
                    %if bound crossing happens 
                choose_strg_tb = nan(max_dur, 2);            
                choose_strg_tb(:,1) = abs(log(Ptb_marginal(:,1,2)./Ptb_marginal(:,1,1))) < theta;     %lower bound crossing
                choose_strg_tb(:,2) = abs(log(Ptb_marginal(:,2,1)./Ptb_marginal(:,2,2))) < theta2;     %upper bound crossing

                    %repeat to set up alternate criterion for stim trials (when using thetaOffset)
                        %if bound crossing does not happen
                choose_strg_xt_stim = nan(length(xmesh), max_dur);
                I = xmesh>0;
                choose_strg_xt_stim(I,:) = abs(log(Pxt_marginal(I,:,1)./Pxt_marginal(I,:,2))) < theta + dtheta;
                I = xmesh<0;
                choose_strg_xt_stim(I,:) = abs(log(Pxt_marginal(I,:,2)./Pxt_marginal(I,:,1))) < theta2 + dtheta2;
                I = xmesh==0;
                if sum(I)==1
                    choose_strg_xt_stim(I,:) = mean([abs(log(Pxt_marginal(find(I)+1,:,2)./Pxt_marginal(find(I)+1,:,1))) ; abs(log(Pxt_marginal(find(I)-1,:,2)./Pxt_marginal(find(I)-1,:,1)))]) < theta;
                end
                    %if bound crossing happens 
                choose_strg_tb_stim = nan(max_dur, 2);            
                choose_strg_tb_stim(:,1) = abs(log(Ptb_marginal(:,1,2)./Ptb_marginal(:,1,1))) < theta + dtheta;     %lower bound crossing
                choose_strg_tb_stim(:,2) = abs(log(Ptb_marginal(:,2,1)./Ptb_marginal(:,2,2))) < theta2 + dtheta2;     %upper bound crossing
        
            %plot the region
        if opt.plot
            plot_logodds(Pxt_marginal, Ptb_marginal, theta, choose_strg_tb, xmesh);
        end
        
            %add a new field to trial_data to make sure that we don't repeat the calculations for
            %trials with similar duration, and motion strength
        trial_data.processed = zeros(size(trial_data.Coh));
        trial_data.expectedPright = nan(size(trial_data.Coh));
        trial_data.expectedPS = nan(size(trial_data.Coh));
        
        try
            n = 0;
            LL = 0;
            t = (1 : max_dur)';
            for s = 0 : 1
                for c = 1 : length(coh_set)
                    Pxt = Pxt_coh(:,:,s+1,c);
                    Ptb = Ptb_coh(:,:,s+1,c);
                    for Ts = 0 : 1
                            %calculate Pright, Pleft and Pstrg
                        if Ts == 0
                            if opt.use_midline_as_criterion==false
                                chance = cumsum(Ptb(:,2).*(LogPostOdds_tb(:,2)==0)) + ...
                                         cumsum(Ptb(:,1).*(LogPostOdds_tb(:,1)==0)) + ...
                                         sum(Pxt.*(LogPostOdds_xt==0),1)';
                                Pright = cumsum(Ptb(:,2).*(LogPostOdds_tb(:,2)>0)) + ...
                                         cumsum(Ptb(:,1).*(LogPostOdds_tb(:,1)>0)) + ...
                                         sum(Pxt.*(LogPostOdds_xt>0),1)' + ...
                                         0.5*chance;
                                Pleft  = cumsum(Ptb(:,2).*(LogPostOdds_tb(:,2)<0)) + ...
                                         cumsum(Ptb(:,1).*(LogPostOdds_tb(:,1)<0)) + ...
                                         sum(Pxt.*(LogPostOdds_xt<0),1)' + ...
                                         0.5*chance;
                            else
                                Pright = cumsum(Ptb(:,2)) + sum(Pxt(xmesh>0,:),1)' + 0.5*Pxt(xmesh==0,:)';
                                Pleft  = cumsum(Ptb(:,1)) + sum(Pxt(xmesh<0,:),1)' + 0.5*Pxt(xmesh==0,:)';
                            end
                            Pstrg = zeros(size(t));
                        elseif Ts == 1
                            if s == 0
                                if opt.use_midline_as_criterion==false
                                    chance = cumsum(Ptb(:,2).*(LogPostOdds_tb(:,2)==0).*(choose_strg_tb(:,2)==0)) + ...
                                             cumsum(Ptb(:,1).*(LogPostOdds_tb(:,1)==0).*(choose_strg_tb(:,1)==0)) + ...
                                             sum(Pxt.*(LogPostOdds_xt==0).*(choose_strg_xt==0))';
                                    Pright = cumsum(Ptb(:,2).*(LogPostOdds_tb(:,2)>0).*(choose_strg_tb(:,2)==0)) + ...
                                             cumsum(Ptb(:,1).*(LogPostOdds_tb(:,1)>0).*(choose_strg_tb(:,1)==0)) + ...
                                             sum(Pxt.*(LogPostOdds_xt>0).*(choose_strg_xt==0))' + ...
                                             0.5*chance;
                                    Pleft  = cumsum(Ptb(:,2).*(LogPostOdds_tb(:,2)<0).*(choose_strg_tb(:,2)==0)) + ...
                                             cumsum(Ptb(:,1).*(LogPostOdds_tb(:,1)<0).*(choose_strg_tb(:,1)==0)) + ...
                                             sum(Pxt.*(LogPostOdds_xt<0).*(choose_strg_xt==0))' + ...
                                             0.5*chance;
                                    Pstrg  = cumsum(Ptb(:,2).*(choose_strg_tb(:,2)==1)) + ...
                                             cumsum(Ptb(:,1).*(choose_strg_tb(:,1)==1)) + ...
                                             sum(Pxt.*(choose_strg_xt==1))';
                                else
                                    Pright = cumsum(Ptb(:,2).*(choose_strg_tb(:,2)==0)) + ...
                                             sum(Pxt(xmesh>0,:).*(choose_strg_xt(xmesh>0,:)==0),1)' + ...
                                             0.5*(Pxt(xmesh==0,:).*(choose_strg_xt(xmesh==0,:)==0))';
                                    Pleft  = cumsum(Ptb(:,1).*(choose_strg_tb(:,1)==0)) + ...
                                             sum(Pxt(xmesh<0,:).*(choose_strg_xt(xmesh<0,:)==0),1)' + ...
                                             0.5*(Pxt(xmesh==0,:).*(choose_strg_xt(xmesh==0,:)==0))';
                                    Pstrg  = cumsum(Ptb(:,2).*(choose_strg_tb(:,2)==1)) + ...
                                             cumsum(Ptb(:,1).*(choose_strg_tb(:,1)==1)) + ...
                                             sum(Pxt.*(choose_strg_xt==1),1)';
                                end
                            else
                                if opt.use_midline_as_criterion==false
                                    chance = cumsum(Ptb(:,2).*(LogPostOdds_tb(:,2)==0).*(choose_strg_tb_stim(:,2)==0)) + ...
                                             cumsum(Ptb(:,1).*(LogPostOdds_tb(:,1)==0).*(choose_strg_tb_stim(:,1)==0)) + ...
                                             sum(Pxt.*(LogPostOdds_xt==0).*(choose_strg_xt_stim==0))';
                                    Pright = cumsum(Ptb(:,2).*(LogPostOdds_tb(:,2)>0).*(choose_strg_tb_stim(:,2)==0)) + ...
                                             cumsum(Ptb(:,1).*(LogPostOdds_tb(:,1)>0).*(choose_strg_tb_stim(:,1)==0)) + ...
                                             sum(Pxt.*(LogPostOdds_xt>0).*(choose_strg_xt_stim==0))' + ...
                                             0.5*chance;
                                    Pleft  = cumsum(Ptb(:,2).*(LogPostOdds_tb(:,2)<0).*(choose_strg_tb_stim(:,2)==0)) + ...
                                             cumsum(Ptb(:,1).*(LogPostOdds_tb(:,1)<0).*(choose_strg_tb_stim(:,1)==0)) + ...
                                             sum(Pxt.*(LogPostOdds_xt<0).*(choose_strg_xt_stim==0))' + ...
                                             0.5*chance;
                                    Pstrg  = cumsum(Ptb(:,2).*(choose_strg_tb_stim(:,2)==1)) + ...
                                             cumsum(Ptb(:,1).*(choose_strg_tb_stim(:,1)==1)) + ...
                                             sum(Pxt.*(choose_strg_xt_stim==1))';
                                else
                                    Pright = cumsum(Ptb(:,2).*(choose_strg_tb_stim(:,2)==0)) + ...
                                             sum(Pxt(xmesh>0,:).*(choose_strg_xt_stim(xmesh>0,:)==0),1)' + ...
                                             0.5*(Pxt(xmesh==0,:).*(choose_strg_xt_stim(xmesh==0,:)==0))';
                                    Pleft  = cumsum(Ptb(:,1).*(choose_strg_tb_stim(:,1)==0)) + ...
                                             sum(Pxt(xmesh<0,:).*(choose_strg_xt_stim(xmesh<0,:)==0),1)' + ...
                                             0.5*(Pxt(xmesh==0,:).*(choose_strg_xt_stim(xmesh==0,:)==0))';
                                    Pstrg  = cumsum(Ptb(:,2).*(choose_strg_tb_stim(:,2)==1)) + ...
                                             cumsum(Ptb(:,1).*(choose_strg_tb_stim(:,1)==1)) + ...
                                             sum(Pxt.*(choose_strg_xt_stim==1),1)';
                                end                                
                            end
                            
                        end
                            %normalize prob to ensure they sum up to one
                        Ptot = Pright + Pleft + Pstrg;
                        if any(abs(Ptot-1)>1e-3)
                            warning('the probabilities (Pright, Pleft, and Pstrg) do not add up to one!\n\tPtot=%f (for s:%d, c:%d, Ts:%d)\n', Ptot, s, c, Ts);
                        end
                        Pright = Pright./Ptot;
                        Pleft = Pleft./Ptot;
                        Pstrg = Pstrg./Ptot;
                            %mark the trials that should be processed now
                        I = trial_data.uStim==s & trial_data.Coh==coh_set(c) & trial_data.SureT==Ts;
                            %calculate the expected probability correct and probability sure target 
                        dur = trial_data.StimDur(I);
                        trial_data.expectedPright(I) = Pright(dur)./(Pright(dur)+Pleft(dur));
                        trial_data.expectedPS(I) = Pstrg(dur);

                            %update log-likelihood 
                        R = unique(trial_data.CorrTrg(trial_data.Coh>0));    %which choice is the rightward choice
                        dur_right = trial_data.StimDur(I & trial_data.Resp==R);
                        dur_left = trial_data.StimDur(I & trial_data.Resp==3-R);
                        dur_strg = trial_data.StimDur(I & trial_data.Resp==3);
                        minP = 1e-300;
                        switch opt.fitwhat{s+1}
                            case 'PrightPS'
                                    % fit Pright and PS
                                LL = LL + sum(log(max(Pstrg(dur_strg),minP))) + sum(log(max(Pright(dur_right),minP))) + sum(log(max(Pleft(dur_left),minP)));
                                n = n + length(dur_right) + length(dur_left) + length(dur_strg);
                            case 'Pright'
                                    % fit only Pright. Renormalize Pright for trials in which Ts was shown
                                Pright = Pright./(Pright+Pleft);
                                LL = LL + sum(log(max(Pright(dur_right),minP))) + sum(log(max(1-Pright(dur_left),minP)));
                                n = n + length(dur_right) + length(dur_left);
                            case 'PS'
                                    % fit Pright for trials without Ts, and only Psure for trials with Ts
                                    % do not fit Pright for trials in which Ts was shown 
                                if Ts == 0
                                    LL = LL + sum(log(max(Pright(dur_right),minP))) + sum(log(max(Pleft(dur_left),minP)));
                                    n = n + length(dur_right) + length(dur_left);
                                elseif Ts == 1
                                    LL = LL + sum(log(max(Pstrg(dur_strg),minP))) + sum(log(max(1-P_strg(dur_right),minP))) + sum(log(max(1-P_strg(dur_left),minP)));
                                    n = n + length(dur_right) + length(dur_left) + length(dur_strg);
                                end
                            case 'Pright-woTs'
                                    %fit only Pright and only for trials in which Ts was not shown
                                if Ts==0
                                        %note that we will not use trials with Ts for this
                                        %fit, multiply LLE by 2 to compensate for the lower
                                        %number of trials in this condition 
                                    LL = LL + 2*(sum(log(max(Pright(dur_right),minP))) + sum(log(max(Pleft(dur_left),minP))));
                                    n = n + length(dur_right) + length(dur_left);
                                end
                            otherwise
                                error('unknown fitwhat!');
                        end
                    end
                end
            end
        catch err
            fprintf('\n\nERROR!!!\n');
            fprintf('run %d\n', call_num);
            fprintf('\tk= %g\n\tB= %g\n\ttheta= %g\n\tdcoh_stim= %g\n\tdcoh_bias= %g\n\tdk= %g\n\tdsigma= %g\n\tweber_fraction= %f\n\t', ...
                    k, B, theta, dcoh_stim, dcoh_bias, dk, dsigma, weber_fraction);
            fprintf('s:%d, c:%d, Ts:%d\n', s, c, Ts);
            for i = 1 : length(err.stack)
                fprintf('entry %d of error stack\n', i);
                fprintf('\tfile: %s\n', err.stack(i).file);
                fprintf('\tfunction: %s\n', err.stack(i).name);
                fprintf('\tline: %d\n', err.stack(i).line);
            end
            rethrow(err)
        end
        
        err = -LL;
        
        % print progress report!
        fprintf ( '\n\n\n****************************************\n');
%         fprintf ( 'file: %s\n' , filename );
        fprintf ( 'run %d\n' , call_num);
        fprintf ( '\tk= %f,\n\tB= %f\n\ttheta= %f,\n\tdcoh_stim= %f\n\tdcoh_bias= %f\n\tdk= %f\n\tdsigma= %f\n\tsigma2beta= %f\n\ttheta2= %f\n\tdtheta= %f\n\tdtheta2= %f\n', k, B, theta, dcoh_stim, dcoh_bias, dk, dsigma, sigma2_beta, theta2, dtheta, dtheta2);
        fprintf ( 'err: %f\n' , err );
        fprintf ('number of processed trials: %d\n', n);
        
%         %print progress report!
%         fprintf('\n\n\n****************************************\n');
%         fprintf('run %d\n', call_num);
%         fprintf('\tk= %g\n\tB= %g\n\ttheta= %g\n\tdcoh_stim= %g\n\tdcoh_bias= %g\n\tdk= %g\n\tdsigma= %g\n\tweber_fraction= %f\n\t', ...
%                 k, B, theta, dcoh_stim, dcoh_bias, dk, dsigma, weber_fraction);
%         fprintf('err: %f\n', err);
%         fprintf('number of processed trials: %d\n', n);
        
        if feedback>0
            if fh>0
                figure(fh);
                plot(call_num, err, '.');
                drawnow;
            end
        end
        call_num = call_num + 1;
    end
end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
            %this function retrieves the full parameter set given the adjustable and
            %fixed parameters
function param2 = getParam ( param1 , guess , fixed )
    param2(fixed==0) = param1;              %get adjustable parameters from param1
    param2(fixed==1) = guess(fixed==1);     %get fixed parameters from guess
end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









