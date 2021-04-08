function processKinarmData(projectPath,combineSameBlocks,plotTrials)
% processKinarmData Process and plot raw Kinarm data trial by trial
%
% processKinarmData loads and processed raw Kinarm data from loadKinarmData.m.
% Raw data needs to be in a folder named 1_RawData_mat, with a folder for
% each subject.
%
% For each trial:
% Get the stimuli and timing of events in the task.
% Get the hand data, low-pass filter acceleration.
% Get the gaze data, remove blinks, low-pass filter, compute angular gaze
%   velocity for the detection of saccades.
% Plot the hand and gaze data (optional).
% Perform custom processing by calling customKinarmTrialProcessing.m.
% Save data in struct E (experiment info), D (data aligned to target
%   appearance) and C (custom variables) in folder 2_ProcessedData.
% 
% processKinarmData(projectPath) allows to define the path with the data of
%   the Kinarm project.
% 
% processKinarmData(projectPath,combineSameBlocks) combines blocks of the
%   same task into a single data file if combineSameBlocks=true 
%   (default=true).
%
% processKinarmData(projectPath,combineSameBlocks,plotTrials) creates a 
%   plot with the hand and gaze data for each trial if plotTrials=true
%   (default=false).
%
% This code was written for Miriam Spering's lab at UBC Vancouver, BC Canada.

% MIT License
% Copyright (c) 2021 Anouk de Brouwer

close all;

% set project path if it isn't defined
if nargin==0 || isempty(projectPath)
    disp('Select folder with KINARM project')
    projectPath = uigetdir;
end

% set defaults
if nargin<=1
    combineSameBlocks = true;
    plotTrials = false;
end

%% Select data files

% check if we are at the right level
while ~exist([projectPath filesep '1_RawData_mat'],'dir');
    expFolder = selectFiles([projectPath '*'],'folders');
    projectPath = [projectPath filesep expFolder.name filesep];
end
% define path for input and output data and figures
dataPath = [projectPath filesep '1_RawData_mat' filesep];
saveToPath = [projectPath filesep '2_ProcessedData' filesep];
saveFigsToPath = [projectPath filesep 'Figures' filesep 'Individual' filesep];
if ~exist(saveToPath,'dir')
    mkdir(saveToPath)
end
if ~exist(saveFigsToPath,'dir')
    mkdir(saveFigsToPath)
end

% select subjects
subjFolders = selectFiles([dataPath '*'],'folders');
nSubj = length(subjFolders);

% load text file with custom trial specs
custom_specs = [];
fileID = fopen([dataPath 'custom_trial_specs.txt'],'r');
if fileID>0
    text = fscanf(fileID,'%s');
    custom_specs = strsplit(text,{',', ', '});
    custom_specs = custom_specs(~cellfun(@isempty,custom_specs));
    fclose(fileID);
end

% load text file with visual event names
events_visual = [];
fileID = fopen([dataPath 'visual_event_names.txt'],'r');
if fileID>0
    text = fscanf(fileID,'%s');
    events_visual = strsplit(text,{',', ', '});
    events_visual = events_visual(~cellfun(@isempty,events_visual));
    fclose(fileID);
end

% open figures
if plotTrials
    fig1 = scaledFigure(2,1);
    set(fig1,'visible','off');
end
fig2 = scaledFigure(2,2);
colors = get(gca,'colororder'); % get default colors
set(fig2,'visible','off');

%% Loop over subjects

for s = 1 : nSubj
    
    % select data files
    dataFiles = selectFiles([dataPath subjFolders(s).name],'files');
    dataFileNames = {dataFiles(:).name}';
    i_ = strfind(dataFileNames{1},'_');
    blockNames = cellfun(@(x) x(i_(2)+1:end-4),dataFileNames,'UniformOutput',false);
    nBlocks = length(dataFiles);
    taskBlockNo = ones(nBlocks,1);
    
    % if blocks of the same task are combined, get number of task block
    if combineSameBlocks
        for b = 1 : nBlocks
            i_ = strfind(dataFileNames{b},'_'); % subjID_expDate_taskProtocol(_blockNumber).mat
            taskBlockNo_str = dataFileNames{b}(i_(end)+1:end-4);
            if isstrprop(taskBlockNo_str,'digit')
                taskBlockNo(b) = str2double(taskBlockNo_str);
            end
        end
    end
    
    %% Loop over blocks
    
    for b = 1 : nBlocks
        
        % create structs for saving of data
        if taskBlockNo(b)==1
            C = [];
            D = [];
            E.calibration = [];
            E.experiment = [];
        end
        
        %% LOAD DATA
        
        % load data file
        data_block = [];
        load([dataPath subjFolders(s).name '/' dataFiles(b).name])
        fprintf('\n')
        disp(['Loading ' dataFiles(b).name '...'])
        nTrials = length(data_block);
        
        % get setup info
        if s==1 && b==1
            fs = expInfo.hand.RATE; % sampling frequency
            displayDelay = 32; % temp
            eyeHeight = 0.35; % temp
            %displayDelay = expInfo.experiment.DISPLAY_DELAY_MS; % display delay (measured using photodiode)
            %eyeHeight = expInfo.experiment.EYE_HEIGHT_M;
        end
        
        % get trial types and number of repetitions from block table
        blockTable = expInfo.block_table.TP_LIST;
        blockTable = blockTable(cellfun(@length,blockTable)>0);
        list = [];
        for r = 1 : length(blockTable)
            str = expInfo.block_table.TP_LIST{r};
            str = strrep(str,'-',':');
            l = eval(['[' str ']']);
            n = expInfo.block_table.LIST_REPS(r)*expInfo.block_table.BLOCK_REPS(r);
            list = [list; repmat(l(:),n,1)];
        end
        iCurrTypes = unique(list);
        nCurrTypes = length(iCurrTypes);
        nRepsType = zeros(1,length(iCurrTypes));
        for t = 1 : nCurrTypes
            nRepsType(t) = sum(list==iCurrTypes(t));
        end
        nReps = max(nRepsType);
        
        % pre-allocate variables
        % repetitions in rows, trial types in columns
        D_block.trialType   = NaN(nReps,nCurrTypes);
        D_block.trialNumber = NaN(nReps,nCurrTypes);
        D_block.tEvents     = cell(nReps,nCurrTypes);
        D_block.xHand       = cell(nReps,nCurrTypes);
        D_block.yHand       = cell(nReps,nCurrTypes);
        D_block.vxHand      = cell(nReps,nCurrTypes);
        D_block.vyHand      = cell(nReps,nCurrTypes);
        D_block.accxHand    = cell(nReps,nCurrTypes);
        D_block.accyHand    = cell(nReps,nCurrTypes);
        D_block.xGaze       = cell(nReps,nCurrTypes);
        D_block.yGaze       = cell(nReps,nCurrTypes);
        D_block.wGaze       = cell(nReps,nCurrTypes);
        D_block.blink       = cell(nReps,nCurrTypes);
        D_block.xGaze_preTarget = cell(nReps,nCurrTypes);
        D_block.yGaze_preTarget = cell(nReps,nCurrTypes);
        D_block.wGaze_preTarget = cell(nReps,nCurrTypes);
        D_block.blink_preTarget = cell(nReps,nCurrTypes);
        C_block = [];
        
        %% BLOCK SPECS - Get trial, stimulus and event specifications
        
        % optional: define custom trial specs from TP table that you want
        % to save in the data file (if they haven't been defined yet)
        while isempty(custom_specs)
            specs = fieldnames(expInfo.tp_table);
            specsNames = [num2cell((1:length(specs))') specs] % show specs names with numbers
            specsNumbers = input(['Custom trial specs are not defined for this task.\n',...
                'Select custom specs by typing the numbers (e.g., [1:5,7]): ']);
            custom_specs = specs(specsNumbers);
            disp('You selected the following custom trial specs:')
            disp(custom_specs)
            saveText = input('Would you like to save this selection for future use? yes(1) or no(0): ');
            % save text file with custom trial specs
            if saveText
                fileID = fopen([dataPath 'custom_trial_specs.txt'],'w');
                fprintf(fileID,'%s,',custom_specs{:});
                fclose(fileID);
            end
        end
        
        % extract trial specifics for each trial type
        trialSpecs = getKinarmTrialSpecs(expInfo,custom_specs);
        
        % get xy range of stimuli
        stimX_min = min([trialSpecs.xStart-trialSpecs.startRadius;...
            trialSpecs.xTarget-trialSpecs.targetRadius;...
            trialSpecs.xTarget-0.5*trialSpecs.targetWidth]);
        stimX_max = max([trialSpecs.xStart+trialSpecs.startRadius;...
            trialSpecs.xTarget+trialSpecs.targetRadius;...
            trialSpecs.xTarget+0.5*trialSpecs.targetWidth]);
        stimY_min = min([trialSpecs.yStart-trialSpecs.startRadius;...
            trialSpecs.yTarget-trialSpecs.targetRadius;...
            trialSpecs.yTarget-0.5*trialSpecs.targetHeight]);
        stimY_max = max([trialSpecs.yStart+trialSpecs.startRadius;...
            trialSpecs.yTarget+trialSpecs.targetRadius;...
            trialSpecs.yTarget+0.5*trialSpecs.targetHeight]);
        stimXY_range = [stimX_min stimX_max stimY_min stimY_max]/100; % in m
        
        % get trial specs for current block only
        trialSpecs = structfun(@(x) x(iCurrTypes,:)',trialSpecs,'UniformOutput',false);
        
        % get all task events
        events = expInfo.event_definitions.LABELS;
        
        % define visual task events if they haven't been defined yet
        while isempty(events_visual)
            eventNames = [num2cell((1:length(events))') events(:)] % show event names with numbers
            eventNumbers = input(['Visual events are not defined for this task.\n',...
                'Select visual events by typing the numbers (e.g., [1:5,7]): ']);
            events_visual = events(eventNumbers);
            disp('You selected the following visual events:')
            disp(events_visual)
            saveText = input('Would you like to save this selection for future use? yes(1) or no(0): ');
            % save text file with visual task events
            if saveText
                fileID = fopen([dataPath 'visual_event_names.txt'],'w');
                fprintf(fileID,'%s,',events_visual{:});
                fclose(fileID);
            end
        end
        eventIsVisual = ismember(events,events_visual);
        
        %% Loop over trials
        
        for t = 1 : nTrials
            
            % get data
            T = data_block(t);
            
            %% TRIAL TYPE - get trial type and number
            
            tp = T.TRIAL.TP;
            tpcol = find(tp==iCurrTypes);
            rep = T.TRIAL.TP_NUM;
            if ~isnan(D_block.trialType(rep,tpcol)) % add rep when block was restarted
                rep = find(isnan(D_block.trialType(:,tpcol)),1);
            end
            D_block.trialType(rep,tpcol) = tp;
            D_block.trialNumber(rep,tpcol) = T.TRIAL.TRIAL_NUM;
            
            %% STIMULI - get position and size of start and target
            % current trial only, in m
            
            % start position and size
            stim.startX = trialSpecs.xStart(tpcol)/100;
            stim.startY = trialSpecs.yStart(tpcol)/100;
            stim.startRadius = trialSpecs.startRadius(tpcol)/100;
            
            % reach target position and size
            stim.targetShape = trialSpecs.targetShape{tpcol};
            stim.targetX = trialSpecs.xTarget(end,tpcol)/100;
            stim.targetY = trialSpecs.yTarget(end,tpcol)/100;
            stim.targetDistance = stim.targetY - stim.startY;
            stim.targetX_h = trialSpecs.xTarget_haptic(end,tpcol)/100; % haptic target position
            stim.targetY_h = trialSpecs.yTarget_haptic(end,tpcol)/100;
            if strcmp(stim.targetShape,'circle')
                stim.targetRadius = trialSpecs.targetRadius(1,tpcol)/100;
                stim.targetWidth = NaN;
                stim.targetHeight = NaN;
            elseif strcmp(stim.targetShape,'rect')
                stim.targetRadius = NaN;
                stim.targetWidth = trialSpecs.targetWidth(end,tpcol)/100;
                stim.targetHeight = trialSpecs.targetHeight(end,tpcol)/100;
            end
            % all target positions
            allTargetXY = unique([trialSpecs.xTarget(end,:)' trialSpecs.yTarget(end,:)'],'rows');
            allTargetXY = allTargetXY(~isnan(allTargetXY(:,1)),:)/100;
            
            % cursor size and visibility
            stim.cursorRadius = trialSpecs.cursorRadius(tpcol);
            stim.cursorVisible = trialSpecs.cursorVisible(tpcol);
            
            %% TRIAL TIMING - get timing of events
            
            % get events in current trial
            events_trial = cellfun(@(x) x(isstrprop(x,'graphic')),...
                T.EVENTS.LABELS,'UniformOutput',false); % trim white space
            
            % check for timeout
            col = getColumnIndex(events_trial,'TIMEOUT');
            if col>0;
                C_block.error{rep,tpcol} = 'Timeout';
                continue
            end
            
            % pre-allocate timing variables
            tEvent = NaN(1,length(events));
            tSend = NaN(1,length(events));
            tDisplay = NaN(1,length(events));
            tTrial = NaN(1,length(events));
            
            % get timing of events
            for e = 1 : length(events) % loop over task events
                col = getColumnIndex(events_trial,events(e)); % get corresponding event in trial
                if col>0
                    % get event time
                    tEvent(e) = T.EVENTS.TIMES(col(end));
                    % get send to display time and video acknowlegdment time for visual events
                    if eventIsVisual(e)
                        i = find(expInfo.video_latency.SEND_TIMES>=tEvent(e),1); % temp
                        %i = find(T.VIDEO_LATENCY.SEND_TIMES>=tEvents(e),1);
                        if ~isempty(i)
                            tSend(e) = expInfo.video_latency.SEND_TIMES(i); % temp
                            tDisplay(e) = expInfo.video_latency.ACK_TIMES(i) + displayDelay/1000; % temp
                            %tSend(e) = T.VIDEO_LATENCY.SEND_TIMES(i);
                            %tDisplay(e) = T.VIDEO_LATENCY.ACK_TIMES(i) + displayDelay/1000; % add display delay
                        end
                    end
                end
            end
            
            % get timing of target onset
            col = getColumnIndex(events,'TARGET_ON'); % modify if event is named differently
            tTarget = tDisplay(col); % use the time that the target appeared on the display
            % temp - if video_latency wasn't saved correctly
            if isnan(tTarget)
                tTarget = tEvent(col) + 57/1000; % event time + total delay
                keyboard
            end
            %
            tTarget_ms = round(tTarget*1000);
            
            % get timing of target offset
            tTargetOff = NaN;
            col = getColumnIndex(events,'TARGET_OFF'); % modify if event is named differently
            if col>0;
                tTargetOff = tDisplay(col); % display time
            end
            % temp - if video_latency wasn't saved correctly
            if isnan(tTargetOff)
                tTargetOff = tEvent(col) + 57/1000;
            end
            
            % define end of trial as target off or end of recording,
            % whichever comes first
            tEnd_ms = min([round(tTargetOff*1000) length(T.Right_HandX)]);
            
            % align timing of visual and non-visual events to target onset
            tTrial(eventIsVisual) = tDisplay(eventIsVisual) - tTarget;
            tTrial(~eventIsVisual) = tEvent(~eventIsVisual) - tTarget;
            tTrial_ms = round(tTrial*1000) + 1; % in ms, target onset at t=1
            
            %% HAND DATA - get position, velocity, and acceleration
            
            % get hand position and remove any 'spikes' in the data
            xyHand_raw = [T.Right_HandX(tTarget_ms:tEnd_ms) T.Right_HandY(tTarget_ms:tEnd_ms)];
            xyHand = medfilt1(xyHand_raw,3); % remove spikes (third-order median filtering)
            d_medfilt = abs(xyHand-xyHand_raw);
            if max(d_medfilt(:)>0.001) % check if any spikes were removed
                plot([xyHand_raw,xyHand]); disp('Removed spike'); keyboard
            end
            
            % get hand velocity
            vxyHand = [T.Right_HandXVel(tTarget_ms:tEnd_ms) T.Right_HandYVel(tTarget_ms:tEnd_ms)];
            
            % get and filter hand acceleration
            [bbutter,abutter] = butter(2,50/(fs/2)); % 2nd order 50 Hz low-pass butterworth filter
            accx = filtfilt(bbutter,abutter,T.Right_HandXAcc);
            accy = filtfilt(bbutter,abutter,T.Right_HandYAcc);
            accxyHand = [accx(tTarget_ms:tEnd_ms) accy(tTarget_ms:tEnd_ms)];
            
            % get estimate of cursor position on screen
            xyCursor_est = xyHand + vxyHand*expInfo.video_settings.FEED_FORWARD;
            
            %% GAZE DATA - remove blinks, low-pass filter, compute angular velocity
            
            % get and filter gaze data
            if expInfo.calibration.GAZE_CALIBRATED == 1
                
                xyLimits = [-25 25 0 50]; % approximate size of the workspace in cm
                nsMax = 100; % maximum number of continuous samples with pupil loss for which gaze position is interpolated
                [xyGaze_filt,vxyGaze_filt,paddedBlink,notUpdated] = filterGazeData(...
                    [T.Gaze_X T.Gaze_Y]*100,xyLimits,fs,nsMax,T.Gaze_TimeStamp);
                if any(notUpdated); disp('Gaze not updated'); keyboard; end % temp
                
                % 2D Cartesian coordinates to 3D eye-centered Cartesian coordinates
                xyzGaze = xyGaze_filt/100; % in m
                xyzGaze(~isnan(xyzGaze(:,1)),3) = eyeHeight;
                vxyzGaze = vxyGaze_filt/100; % in m/s
                vxyzGaze (~isnan(vxyzGaze(:,1)),3) = 0; % assume eye is stationary in z
                
                % 3D Cartesian coordinates to spherical coordinates, and
                % spherical coordinates to angular gaze velocity -
                % angular velocity should be used for saccade detection
                vGaze_ang = cart2angularVelocity(xyzGaze,vxyzGaze);
                
                % select data from target onset until end of trial
                xyGaze = xyGaze_filt(tTarget_ms:tEnd_ms,:)/100; % in m
                vGaze_angular = vGaze_ang(tTarget_ms:tEnd_ms); % in deg/s
                blink = paddedBlink(tTarget_ms:tEnd_ms);
                
                % select data before target onset
                % this can be useful to detect early/anticipatory saccades
                n = 500;
                tPre_ms = max([tTarget_ms-n 1]);
                int_pre = tPre_ms:tTarget_ms-1;
                xyGaze_preTarget = NaN(n,2);
                vGaze_angular_preTarget = NaN(n,1);
                blink_preTarget = NaN(n,1);
                xyGaze_preTarget(int_pre-max(int_pre)+n,:) = xyGaze_filt(int_pre,:)/100; % in m;
                vGaze_angular_preTarget(int_pre-max(int_pre)+n,:) = vGaze_ang(int_pre,:); % in deg/s
                blink_preTarget(int_pre-max(int_pre)+n,:) = paddedBlink(int_pre,:);
                
            else % no gaze data
                xyGaze = [NaN NaN];
                vGaze_angular = NaN;
                blink = NaN;
                xyGaze_preTarget = [NaN NaN];
                vGaze_angular_preTarget = NaN;
                blink_preTarget = NaN;
            end
            
            %% STORE - store data of current trial in struct D_block
            
            D_block.tEvents{rep,tpcol} = tTrial_ms(:);
            D_block.xHand{rep,tpcol} = xyHand(:,1);
            D_block.yHand{rep,tpcol} = xyHand(:,2);
            D_block.vxHand{rep,tpcol} = vxyHand(:,1);
            D_block.vyHand{rep,tpcol} = vxyHand(:,2);
            D_block.accxHand{rep,tpcol} = accxyHand(:,1);
            D_block.accyHand{rep,tpcol} = accxyHand(:,2);
            D_block.xGaze{rep,tpcol} = xyGaze(:,1);
            D_block.yGaze{rep,tpcol} = xyGaze(:,2);
            D_block.wGaze{rep,tpcol} = vGaze_angular(:,1);
            D_block.blink{rep,tpcol} = blink;
            D_block.xGaze_preTarget{rep,tpcol} = xyGaze_preTarget(:,1);
            D_block.yGaze_preTarget{rep,tpcol} = xyGaze_preTarget(:,2);
            D_block.wGaze_preTarget{rep,tpcol} = vGaze_angular_preTarget(:,1);
            D_block.blink_preTarget{rep,tpcol} = blink_preTarget;
            
            %% PLOT - plot hand and eye position and velocity
            
            % plot data of current trial
            if plotTrials
                plotKinarmTrial(fig1,rep,tpcol,stim,stimXY_range,...
                    xyHand,xyGaze,vxyHand,vGaze_angular)
            end
            
            %% CUSTOM - custom trial processing defined by the user
            
            % do any custom processing of current trial data and save in C_block
            % e.g., get reach and saccade onsets and offsets and plot them
            % in the current figure, create error variable
            C_block = customKinarmTrialProcessing(C_block,rep,tpcol,plotTrials,expInfo,stim,...
                tTrial_ms,xyHand,vxyHand,accxyHand,xyGaze,vGaze_angular,blink,...
                xyGaze_preTarget,vGaze_angular_preTarget,blink_preTarget);
            
        end % loop over trials
        
        %% ERROR TRIALS - display number of error trials if error variable exists
        
        if isfield(C_block,'error')
            
            % replace empty cells with an empty string
            empty = cellfun(@isempty,C_block.error);
            C_block.error(empty) = {''};
            
            % get number of errors for each type
            uniError = unique(C_block.error);
            uniError = uniError(~(strcmp(uniError,'None') | strcmp(uniError,'')));
            nErrors = zeros(1,length(uniError));
            for e = 1 : length(uniError)
                nErrors(e) = sum(strcmp(C_block.error(:),uniError(e)));
            end
            
            % get number of trials without error
            nGood = sum(strcmp(C_block.error(:),'None'));
            
            % display error table
            varNames = strrep(uniError',' ','_');
            T = array2table([nGood nErrors],'VariableNames',[{'Good'} varNames],...
                'RowNames',{'N_trials'});
            disp(T)
            
            % separate error trials and good trials for plot
            out = ~(strcmp('None',C_block.error) | cellfun(@isempty,C_block.error));
            D_error = replaceErrorTrialsInStruct(D_block,~out,true);
            D_good = replaceErrorTrialsInStruct(D_block,out,true);
            
        end
        
        %% PLOT - plot hand trajectories for each trial type
        
        figName = [subjFolders(s).name '_' blockNames{b}];
        
        figure(fig2); clf
        for c = 1 : nCurrTypes
            subplot(ceil(sqrt(nCurrTypes)),ceil(sqrt(nCurrTypes)),c);
            hold on
            % plot start position
            startXY = [trialSpecs.xStart(c) trialSpecs.yStart(c)];
            startR = [trialSpecs.targetRadius(c) trialSpecs.targetRadius(c)];
            rectangle('Position',[startXY-startR startR*2]/100,'Curvature',[1 1]);
            % plot target position
            if strcmp(trialSpecs.targetShape(c),'circle')
                targetR = [trialSpecs.targetRadius(c) trialSpecs.targetRadius(c)];
                curv = [1 1];
            else
                targetR = 0.5*[trialSpecs.targetWidth(c) trialSpecs.targetHeight(c)];
                curv = [0 0];
            end
            targetXY = [trialSpecs.xTarget(c) trialSpecs.yTarget(c)];
            rectangle('Position',[targetXY-targetR targetR*2]/100,'Curvature',curv)
            % plot hand trajectories
            if exist('D_error','var')
                % plot error trials
                cellfun(@(x,y) plot(x,y,'Color','r'),...
                    D_error.xHand(:,c),D_error.yHand(:,c));
                % plot valid trials
                cellfun(@(x,y) plot(x,y,'Color',colors(1,:),'LineWidth',1),...
                    D_good.xHand(:,c),D_good.yHand(:,c));
            else
                % plot all trials
                cellfun(@(x,y) plot(x,y,'Color',colors(1,:)),...
                    D_block.xHand(:,c),D_block.yHand(:,c));
            end
            % axes and title
            axis equal
            set(gca,'xtick',-0.5:0.1:0.5); set(gca,'ytick',0:0.1:0.5);
            title(['Type ' num2str(D_block.trialType(1,c))])
        end
        suplabel('X (m)','x');
        suplabel('Y (m)','y');
        suplabel(figName,'t');
        saveFigAsPDF([saveFigsToPath 'handTrajectories_' figName],12)
        
        %% SAVE - save block data in struct
        
        % save info about experiment to E struct
        if taskBlockNo(b)==1 % save once
            E.subj = subjFolders(s).name;
            E.trialSpecs = trialSpecs;
            E.accessories = expInfo.accessories;
            E.eventDefinitions = expInfo.event_definitions;
            E.manufacturer = expInfo.manufacturer;
            E.taskWideParams = expInfo.task_wide_params;
            E.videoSettings = expInfo.video_settings;
            E.calibration = expInfo.calibration;
            E.experiment = expInfo.experiment;
        end
        if taskBlockNo>1 % concatenate info for blocks of the same task
            fields = {'calibration','experiment'};
            for f = 1 : length(fields)
                E.(fields{f}) = [E.(fields{f}); expInfo.(fields{f})];
            end
        end
        
        % save data to D struct
        D.samplingFrequency = fs;
        if taskBlockNo(b)==1
            D = D_block;
        else % concatenate blocks of the same task
            fields = fieldnames(D_block);
            fields = fields(~(strcmp(fields,'samplingFrequency') | strcmp(fields,'description')));
            for f = 1 : length(fields)
                D.(fields{f}) = [D.(fields{f}); D_block.(fields{f})];
            end
        end
        % add description
        str = sprintf('All data are sampled at %d Hz and aligned to target onset (t=1), correcting for a display delay (send time to appearance) of %d ms.',...
            fs,displayDelay);
        D.description = ...
            {str;...
            'Timing of events (as defined in expInfo.event_definitions.LABELS) relative to target onset';...
            'Left-right position of handle in m, with screen midline at x=0';...
            'Fore-aft position of handle in m, with y=0 near edge of screen closest to body';...
            'Left-right velocity of handle in m/s';...
            'Fore-aft velocity of handle in m/s';...
            'Left-right acceleration of handle in m/s^2, 50 Hz low-pass filtered';...
            'Fore-aft acceleration of handle in m/s^2, 50 Hz low-pass filtered';...
            'Left-right position of gaze on screen in m';...
            'Fore-aft position of gaze on screen in m';...
            'Angular gaze velocity in deg/s';...
            'Left-right position of gaze on screen, before target onset';...
            'Fore-aft position of gaze on screen, before target onset';...
            'Angular gaze velocity, before target onset'};
        
        % save custom data to C struct
        if taskBlockNo(b)==1
            C = C_block;
        else % concatenate blocks of the same task
            fields = fieldnames(C_block);
            for f = 1 : length(fields)
                C.(fields{f}) = [C.(fields{f}); C_block.(fields{f})];
            end
        end
        
        % save file
        if ~combineSameBlocks || (b==nBlocks || taskBlockNo(b+1)==1)
            fileName = [subjFolders(s).name '_' expInfo.experiment.TASK_PROTOCOL '.mat'];
            if exist([saveToPath fileName],'file') == 2
                disp(['A file named ' fileName ' already exists in ' saveToPath '.'])
                saveFile = input('Do you want to overwrite it? Yes(1) or no(0): ');
                saveFile = logical(saveFile);
            else
                saveFile = true;
            end
            if saveFile == 1
                save([saveToPath fileName],'E','D','C');
                disp(['Saved ' saveToPath fileName])
            else
                disp('Data was not saved')
            end
        end
        
        disp('Press continue to go to next block')
        keyboard
        
    end % loop over blocks
    
    fprintf('\n')
    disp('Press continue to go to next subject')
    keyboard
    
end % loop over subjects

end % end of function


function vGaze_ang = cart2angularVelocity(xyzGaze,vxyzGaze)
% cart2angularVelocity Convert position and velocity in 3D cartesian
% coordinates to angular velocity
%
% for details see: Singh, Perry, Herter (2016). A geometric method for
% computing ocular kinematics and classifying gaze events using monocular
% remote eye-tracking in a robotic environment. J Neuroeng Rehabil 13:10,
% doi 10.1186/s12984-015-0107-4.

% 3D Cartesian coordinates to spherical coordinates
% theta = azimuth ('compass' angle), phi = elevation, in degrees
x = xyzGaze(:,1); y = xyzGaze(:,2); z = xyzGaze(:,3);
[theta,phi,r] = cart2sphd(x,y,z,180);

% compute angular velocity
r_numer = x.*vxyzGaze(:,1) + y.*vxyzGaze(:,2) + 0;
r_dot = r_numer ./ r;
theta_numer = vxyzGaze(:,1).*y - x.*vxyzGaze(:,2);
theta_denom = x.^2 + y.^2;
theta_dot = theta_numer ./ theta_denom; % eq 4b
phi_numer = z.*x.*vxyzGaze(:,1) + z.*y.*vxyzGaze(:,2) - 0;
k = sqrt(x.^2+y.^2);
phi_denom = x.^2.*k + y.^2.*k + z.^2.*k;
phi_dot =  phi_numer ./ phi_denom; % eq 4a
vGaze_ang = sqrt( (theta_dot.*sind(phi)).^2 + phi_dot.^2) / pi*180; % eq 6a, deg/s

end