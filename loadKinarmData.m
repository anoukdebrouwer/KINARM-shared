function loadKinarmData(saveToPath,sortTrialsByExecutionOrder)
% loadKinarmData  First step in processing of raw Kinarm data
%
% loadKinarmData(saveToPath) preprocesses the raw Kinarm data with the
% following steps:
% 1) Load the raw data from the Kinarm using zip_load.m or c3d_load.m (BKIN).
% 2) Get setup and task information (that is the same for each trial).
% 3) Get temporal data and add hand kinematics using KINARM_add_hand_kinematics.m (BKIN).
% 4) Discard fields without data (e.g. left hand for unimanual setup,
%       gaze data, force plate data).
% 5) Save preprocessed data as .mat file in saveToPath/subjID folder.
%       Variable expInfo contains setup and task info, variable data_block
%       contains temporal data.
% Preprocessing can be performed on the Kinarm computer. The preprocessed
% data file does not contain identifying information about the participant.
%
% loadKinarmData(saveToPath,sortTrialsByExecutionOrder) additionally sorts
% trials by execution order if sortTrialsByExecutionOrder is set to true.
% Default is to sort by trial type.

% MIT License
% Copyright (c) 2021 Anouk de Brouwer

% select subject folders for processing
%rawDataPath = 'C:/Users/Bkin/Documents/BKIN Dexterit-E Data/'; % raw data folder on Kinarm computer
rawDataPath = '~/Dropbox/MTLPatients_SHARED/VMRReachPath/'; % temp
cd(rawDataPath)
subj = selectFiles([rawDataPath '*'],'folders');

% select path to save data to if undefined
if nargin==0
    disp('Please select the folder to save the data to')
    saveToPath = uigetdir;
end

% set sortTrialsByExecutionOrder and mergeAfterRestart to false if unspecified
if nargin<=1
    sortTrialsByExecutionOrder = false;
end

% loop over subjects
for s = 1 : length(subj)
    
    %% 1) Load data
    
    % go to subject raw data folder and load data using zip_load
    path_subjRawData = [rawDataPath subj(s).name '/'];
    cd(path_subjRawData)
    fprintf('\nLoading %sRawData/%s...\n',rawDataPath,subj(s).name)
    files = selectFiles(path_subjRawData,'*.zip');
    nBlocks = length(files);
    taskNames = cell(size(files));
    for b = 1 : nBlocks
        D = zip_load(files(b).name); % opens ZIP_FILENAME and outputs the data into the C3D_DATA structure
        c3d_all(b) = D;
        %taskNames{b} = c3d_all(b).c3d(1).EXPERIMENT.TASK_PROGRAM;
        taskNames{b} = c3d_all(b).c3d(1).EXPERIMENT.TASK_PROTOCOL;
    end
    
    % get the number of blocks for each task
    uniTaskNames = unique(taskNames);
    taskBlockNo = zeros(size(taskNames));
    for t = 1 : length(uniTaskNames)
        currTask = strcmp(uniTaskNames{t},taskNames);
        if sum(currTask)>1
            taskBlockNo(currTask) = 1 : sum(currTask);
        end
    end
    
    % loop over blocks
    for b = 1 : nBlocks
        
        % data of current block
        if sortTrialsByExecutionOrder
            c3d_block = sort_trials(c3d_all(b),'execution');
            c3d_block = c3d_block.c3d;
        else
            c3d_block = sort_trials(c3d_all(b),'tp');
            c3d_block = c3d_block.c3d;
        end
        filename = c3d_all(b).file_name;
        
        %% 2) Get setup and task info
        
        % extract information that is the SAME for each trial
        % (non-temporal data fields, all struct fields except EVENTS,
        % TRIAL, and VIDEO_LATENCY)
        c3d_trial1 = c3d_block(1);
        expInfo.hand            = c3d_trial1.HAND;
        expInfo.event_definitions = c3d_trial1.EVENT_DEFINITIONS;
        expInfo.calibration     = c3d_trial1.CALIBRATION;
        expInfo.experiment      = c3d_trial1.EXPERIMENT;
        expInfo.left_kinarm     = c3d_trial1.LEFT_KINARM;
        if isfield(c3d_trial1,'LOAD_TABLE')
            expInfo.load_table  = c3d_trial1.LOAD_TABLE;
        end
        expInfo.task_wide_params = c3d_trial1.TASK_WIDE_PARAMS;
        expInfo.block_table     = c3d_trial1.BLOCK_TABLE;
        expInfo.manufacturer    = c3d_trial1.MANUFACTURER;
        expInfo.right_kinarm    = c3d_trial1.RIGHT_KINARM;
        expInfo.target_table    = c3d_trial1.TARGET_TABLE;
        expInfo.torque_motor    = c3d_trial1.TORQUE_MOTOR;
        expInfo.tp_table        = c3d_trial1.TP_TABLE;
        expInfo.video_settings  = c3d_trial1.VIDEO_SETTINGS;
        expInfo.accessories     = c3d_trial1.ACCESSORIES;
        
        % add info specific to UBC Kinarm setup
        if strcmp(expInfo.experiment.LOCATION,'UBC - Pai')
            expInfo.experiment.EYE_HEIGHT_M = 0.35; % 0.15 m eye to mirror + 0.20 m mirror to screen
            expInfo.experiment.DISPLAY_DELAY_MS = 32; % display delay from photodiode testing
            disp(['Display delay is set to ' num2str(expInfo.experiment.DISPLAY_DELAY_MS) ' ms'])
        end
        
        %% 3) Get temporal data and add hand kinematics
        
        % calculate the hand velocities, accelerations and commanded forces from the
        % joint velocities, accelerations and motor torques for the Kinarm robot
        data_block = KINARM_add_hand_kinematics(c3d_block);
        
        % remove non-temporal fields from data_block structure
        data_block = rmfield(data_block,{'HAND','EVENT_DEFINITIONS',...
            'CALIBRATION','EXPERIMENT','LEFT_KINARM','RIGHT_KINARM',...
            'TASK_WIDE_PARAMS','BLOCK_TABLE','USER_CHANNELS','MANUFACTURER',...
            'TARGET_TABLE','TORQUE_MOTOR','TP_TABLE','VIDEO_SETTINGS','ACCESSORIES'});
        
        %% 4) Discard fields without data
        
        % remove left-hand data if testing on unimanual setup
        if expInfo.left_kinarm.IS_PRESENT == 0
            fieldNames = fieldnames(data_block);
            left = ~cellfun(@isempty,regexp(fieldNames,'Left'));
            data_block = rmfield(data_block,fieldNames(left));
        end
        
        % remove gaze data if eye tracker is not calibrated
        if expInfo.calibration.GAZE_CALIBRATED == 0
            fieldNames = fieldnames(data_block);
            gaze = ~cellfun(@isempty,regexp(fieldNames,'Gaze'));
            data_block = rmfield(data_block,fieldNames(gaze));
        end
        
        % remove force plate data if testing on setup without force plates
        if expInfo.accessories.FORCE_PLATE_COUNT == 0;
            fieldNames = fieldnames(data_block);
            fp = ~cellfun(@isempty,regexp(fieldNames,'FP'));
            data_block = rmfield(data_block,fieldNames(fp));
        end
        
        %% 5) Save data
        
        % create folder for preprocessed data using subject ID
        i_ = strfind(subj(s).name,'_');
        subjID = subj(s).name(i_(end)+1:end); % from folder name
        saveToPath_subj = [saveToPath '/' subjID '/'];
        if ~exist(saveToPath_subj,'dir')
            mkdir(saveToPath_subj)
        end
        
        % create filename: subjID_expDate_taskProtocol(_blockNumber)
        % add block number if there is more than 1 block for this task
        cd(saveToPath_subj)
        %expdate = c3d_trial1.TRIAL.DATE; % from computer date
        expDate = regexp(filename,'2\d+-\d+-\d+','match','once'); % from filename
        filename = [subjID '_' expDate(isstrprop(expDate,'digit')) '_' taskNames{b}];
        if taskBlockNo(b)>0
            filename =  [filename '_' num2str(taskBlockNo(b))];
        end
        
        % check if file does not exist yet
        overwrite = 0;
        while exist([saveToPath_subj filename '.mat'],'file')==2 && overwrite==0
            disp(['A file named ' filename ' already exists.'])
            overwrite = input('Do you want to overwrite it? Yes(1) or no(0): ');
            if overwrite == 0
                filename = input('Provide new filename or press q to return to keyboard: ','s');
                if strcmp(filename,'q')
                    filename = [subjID '_X'];
                    keyboard
                end
            end
        end
        
        % save .MAT file with preprocessed data
        save(filename,'expInfo','data_block')
        disp(['Saved  ' filename])
        
    end % loop over blocks
    
end % loop over subjects