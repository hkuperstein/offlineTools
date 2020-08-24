%% PLDAPS (PDS) cleanup
% after downloading from server, erase all the garbage we don't use.
% this will save disk space, subsequent loading time, and also make
% exploring the data easier.

% VERSION 2.0: 08-30-19 CF
% must be run by getDataFromServer


% important toggles; most or all of these disk-space hogging data 
% can be removed for simple behavioral analyses
removeAnalogData = 1;
removeTimingData = 1;
if strcmp(subject,'human')
    removeMotionTrackingData = 1; % will change this later
    removeDotPositionData = 1;
else
    removeMotionTrackingData = 1;
    removeDotPositionData = 0;
end


fprintf(['\ncleaning up ' remoteFiles{n} '...']);
try
    load([localDir remoteFiles{n}],'-mat');
catch me
    fprintf(' Could not load! Check file. Skipping...\n');
    return
end

try PDS = rmfield(PDS,'initialParameters'); catch; end
try PDS = rmfield(PDS,'initialParameterNames'); catch; end
try PDS = rmfield(PDS,'initialParametersMerged'); catch; end
try PDS = rmfield(PDS,'functionHandles'); catch; end
try PDS = rmfield(PDS,'conditionNames'); catch; end

for t = 1:length(PDS.data) % loop over trials for this file

    if removeAnalogData % also removes analog-derived vars
        try PDS.data{t}.datapixx = rmfield(PDS.data{t}.datapixx,'adc'); catch; end
        try PDS.data{t}.behavior = rmfield(PDS.data{t}.behavior,'fixFP'); catch; end
        try PDS.data{t}.behavior = rmfield(PDS.data{t}.behavior,'eyeXYs'); catch; end
        try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'eyeXYs'); catch; end
    else
        try PDS.data{t}.behavior.eyeXYs = PDS.data{t}.stimulus.eyeXYs; catch; end
        try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'eyeXYs'); catch; end
    end
    
    if removeTimingData
        try PDS.data{t} = rmfield(PDS.data{t},'timing'); catch; end
    end

    if removeMotionTrackingData
        try PDS.data{t} = rmfield(PDS.data{t},'imu'); catch; end
        try PDS.data{t} = rmfield(PDS.data{t},'nexonar'); catch; end
        try PDS.data{t} = rmfield(PDS.data{t},'mp'); catch; end
    end
    
    if removeDotPositionData
        try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'dotX'); catch; end
        try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'dotY'); catch; end
        try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'dotZ'); catch; end
        try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'dotPos'); catch; end        
        try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'dotSize'); catch; end
    end
    
    try PDS.data{t}.behavior = rmfield(PDS.data{t}.behavior,'reward'); catch; end
    try PDS.data{t} = rmfield(PDS.data{t},'mouse'); catch; end
    try PDS.data{t}.behavior.goodtrial = PDS.data{t}.pldaps.goodtrial; catch;  end
    try PDS.data{t} = rmfield(PDS.data{t},'pldaps'); catch; end                               
    try PDS.data{t} = rmfield(PDS.data{t},'plexon'); catch; end
    
    try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'timeFpOn'); catch; end
    try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'timeFpOff'); catch; end
    try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'timeStimOn'); catch; end
    try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'timeStimOff'); catch; end
    try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'timeChoice'); catch; end
    try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'frameFpOn'); catch; end
    try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'frameFpOff'); catch; end
    try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'frameTargetOn'); catch; end
    try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'frameTargetOff'); catch; end
    try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'frameFpOn'); catch; end
    try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'frameFpOff'); catch; end
    try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'frameStimOn'); catch; end
    try PDS.data{t}.stimulus = rmfield(PDS.data{t}.stimulus,'frameStimOff'); catch; end
    
    try PDS.data{t}.behavior.timeTargDisappears = PDS.data{t}.task.timeTargDisappears; catch; end
    try PDS.data{t}.behavior.probOfMemorySaccade = PDS.data{t}.task.probOfMemorySaccade; catch; end
    try PDS.data{t} = rmfield(PDS.data{t},'task'); catch; end
    try PDS.data{t} = rmfield(PDS.data{t},'postTarget'); catch; end                               
    try PDS.data{t} = rmfield(PDS.data{t},'flagNextTrial'); catch; end                               
    try PDS.data{t} = rmfield(PDS.data{t},'iFrame'); catch; end                               
    try PDS.data{t} = rmfield(PDS.data{t},'state'); catch; end                               
    try PDS.data{t} = rmfield(PDS.data{t},'keyboard'); catch; end                               
    try PDS.data{t} = rmfield(PDS.data{t},'currentFrameState'); catch; end                               
    try PDS.data{t} = rmfield(PDS.data{t},'eyeX'); catch; end                               
    try PDS.data{t} = rmfield(PDS.data{t},'eyeY'); catch; end                               
    try PDS.data{t} = rmfield(PDS.data{t},'pupil'); catch; end                               
    try PDS.data{t} = rmfield(PDS.data{t},'remainingFrameTime'); catch; end                               
    try PDS.data{t} = rmfield(PDS.data{t},'framePreLastDrawIdleCount'); catch; end                               
    try PDS.data{t} = rmfield(PDS.data{t},'framePostLastDrawIdleCount'); catch; end                               

    % rename fields (dots3DMP)
    try PDS.conditions{t}.stimulus = renameStructField(PDS.conditions{t}.stimulus,'conflictAngle','delta'); catch; end
    try PDS.conditions{t}.stimulus = renameStructField(PDS.conditions{t}.stimulus,'stimCondition','modality'); catch; end
    try
        if PDS.conditions{t}.stimulus.modality==1 % don't need a bunch of nans for vestib trials, make it just one nan
            PDS.data{t}.stimulus.dotX = NaN;
            PDS.data{t}.stimulus.dotY = NaN;
            PDS.data{t}.stimulus.dotZ = NaN;
            PDS.data{t}.stimulus.dotSize = NaN;
        end
    catch    
    end

end

save([localDir remoteFiles{n}],'PDS','-v7.3');
fprintf(' done\n');

