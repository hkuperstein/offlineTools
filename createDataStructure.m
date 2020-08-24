%% createDataStructure:
% load local files and create a single struct with all the trials and
% fields, excluding fieldExcludes (check struct data after using and add
% excludes as needed)

% VERSION 2.0: 08-30-19 CF
% must be run by PLDAPS_preprocessing


data = struct;
data.filename = {};
data.subj = {};
data.choice = []; % initialize this one field, you'll see why

fieldExcludes = {'reward','leftEarly','fixFP','eyeXYs','corrLoopActive','goodtrial', ...
                 'timeTargDisappears','probOfMemorySaccade','leftTargR','leftTargTheta', ...
                 'rightTargR','rightTargTheta','audioFeedback','textFeedback'};

% now search localDir again for matching files and extract the desired variables from PDS
allFiles = dir(localDir);
for d = 1:length(dateRange)
    for f = 3:length(allFiles) % skip 1+2, they are are "." and ".."
        if contains(allFiles(f).name, subject) ... % check if the target file is in localDir
           && contains(allFiles(f).name, num2str(dateRange(d))) ...
           && contains(allFiles(f).name, paradigm)
       
            disp(['loading ' allFiles(f).name]);
            try
                load([localDir allFiles(f).name],'-mat'); % load it. this is the time limiting step;
                                                          % will eventually change how data are saved to make this faster
                if exist('PDS','var')
                    T = length(data.choice); % set trial counter to continue where left off (0, to start)
                    fprintf('\ncumulative trials processed = %d\n',T);
                    for t = 1:length(PDS.data) % loop over trials for this file,
                        if isfield(PDS.data{t}.behavior,'choice') % and save out the data, excluding trials with
                                                                  % missing data (choice is a good marker for this)
                            T = T+1; % increment trial counter

                            data.filename{T,1} = allFiles(f).name(1:end-4);
                            dateStart = strfind(allFiles(f).name,'20');
                            if contains(subject,'human')
                                data.subj{T,1} = allFiles(f).name(dateStart(1)-3:dateStart(1)-1); % 3-letter code
                            else
                                data.subj{T,1} = subject;
                            end

                            % independent variables are stored in PDS.conditions.stimulus
                            fnames = fieldnames(PDS.conditions{t}.stimulus);
                            fnames(ismember(fnames,fieldExcludes)) = [];
                            for F = 1:length(fnames)
                                eval(['data.' fnames{F} '(T,1) = PDS.conditions{t}.stimulus.' fnames{F} ';']);
                            end
                                % EXCEPT! dot duration, and anything else in
                                % PDS.data.stimulus, i.e. things that are generated
                                % on each trial and not pre-configured. For now
                                % will need to hard-code them here.
                            if isfield(PDS.data{t}.stimulus,'dotDuration')
                                data.duration(T,1) = PDS.data{t}.stimulus.dotDuration;
                            end
                            if isfield(PDS.data{t}.stimulus,'dotPos')
                                data.dotPos{T,1} = PDS.data{t}.stimulus.dotPos;
                            end
                            
                            % dependent variables are stored in PDS.data.behavior
                            fnames = fieldnames(PDS.data{t}.behavior);
                            fnames(ismember(fnames,fieldExcludes)) = [];
                            for F = 1:length(fnames)
                                % SJ 07/20, correct was defaulting to
                                % logical, but then throwing error if a NaN
                                % came up
                                if strcmp(fnames{F},'correct'), data.correct(T,1) = 0; end
                                eval(['data.' fnames{F} '(T,1) = PDS.data{t}.behavior.' fnames{F} ';']);
                            end

                        end
                    end
                    clear PDS
                end
            
            catch
                warning(['error loading ' allFiles(f).name ' -- skipping']);
            end

        end
    end
end


%% some bookkeeping: manually merge any vars whose names have changed, etc

if strcmp(subject(1:5),'human')
    data.conf = data.saccEndPoint;
else
    % do we still need this? SJ 07/2020
    if isfield(data,'postDecisionConfidence')
    if isfield(data,'PDW')
        if length(data.postDecisionConfidence)<length(data.PDW) && length(data.PDW)==length(data.choice)
            data.PDW(1:length(data.postDecisionConfidence)) = data.postDecisionConfidence;
            data = rmfield(data,'postDecisionConfidence');
        else
            error('unsure, diagnose by looking at data');
        end
    elseif length(data.postDecisionConfidence)==length(data.choice)
        data.PDW = data.postDecisionConfidence;
        data = rmfield(data,'postDecisionConfidence');
    else
        error('unsure, diagnose by looking at data');
    end
    end
end
data = rmfield(data,'saccEndPoint'); % either way, this gets removed


disp('done.');

