function trialSpecs = getKinarmTrialSpecs(expInfo,additionalSpecs)
% getKinarmTrialSpecs  Get specifics of trial types from the Kinarm TP table
% and target table
%
% trialSpecs = getKinarmTrialSpecs(expInfo) gets the specifics of all trial
% types from 'expInfo' struct created with loadKinarmData.m and saves it in
% 'trialSpecs' struct. trialSpecs contains arrays with:
% - trial type (type)
% - hand cursor radius (cursorRadius) and visibility (cursorVisible)
% - start position (xStart, yStart) and radius (startRadius)
% - reach target position (xTarget, yTarget), size (targetRadius or
%   targetWidth, targetHeight) and shape (targetShape: 'rect' or 'circle)
%
% trialSpecs = getKinarmTrialSpecs(expInfo,additionalSpecs) saves
% additional specifics that are relevant to the user's experiment.
% Input additionalSpecs should be a cell with the names of the relevant
% variables in the TP table.

% MIT License
% Copyright (c) 2021 Anouk de Brouwer

if nargin==1
    additionalSpecs = [];
end

targetTable = expInfo.target_table;

%% Get TP table and convert to matrix

allfields = fieldnames(expInfo.tp_table);
iFields = find(structfun(@length,expInfo.tp_table)==100);
for i = 1 : length(iFields)
    TPcolnames{i} = allfields{iFields(i)};
    TPtable(:,i) = expInfo.tp_table.(allfields{iFields(i)});
end
notZero = sum(TPtable>0,2)>1; % find rows for which target position is defined
nTP = find(notZero,1,'last');
TPtable = TPtable(notZero,:);
trialSpecs.type(:,1) = 1:nTP;

%% Cursor size and visibility

% pre-allocate
trialSpecs.cursorRadius = NaN(nTP,1);
trialSpecs.cursorVisible = NaN(nTP,1);

% get info
trialSpecs.cursorRadius(notZero,1) = 0.5; % default KINARM value
feedbackMode = expInfo.experiment.HAND_FEEDBACK_MODE;
if strcmp(feedbackMode,'None')
    trialSpecs.cursorVisible(notZero,1) = 0;
elseif ~isempty(strfind(feedbackMode,'hand'))
    trialSpecs.cursorVisible(notZero,1) = 1;
elseif strcmp(feedbackMode,'Controlled by task program')
    % get column index
    iCol = getColumnIndex(TPcolnames,'Cursor_visible'); % modify if variable is named differently
    if iCol>0
        % get info from TP table
        trialSpecs.cursorVisible(notZero,1) = TPtable(:,iCol);
    else
        disp(['Could not find variable that specifies cursor visibility.\n'
            'Make sure the variable is named correctly.'])
        keyboard
    end
end

%% Start target position and size

% get column index
iCol = getColumnIndex(TPcolnames,'Start_target'); % modify if target is named differently

% pre-allocate
trialSpecs.xStart = NaN(nTP,1);
trialSpecs.yStart = NaN(nTP,1);
trialSpecs.startRadius = NaN(nTP,1);

% get info from target table
if iCol>0
    targetNumber = TPtable(:,iCol);
    trialSpecs.xStart(notZero,1) = targetTable.X_GLOBAL(targetNumber);
    trialSpecs.yStart(notZero,1) = targetTable.Y_GLOBAL(targetNumber);
    if isfield(targetTable,'Radius')
        trialSpecs.startRadius(notZero,1) = targetTable.Radius(targetNumber);
    elseif isfield(targetTable,'Visual_radius')
        trialSpecs.startRadius(notZero,1) = targetTable.Visual_radius(targetNumber);
    else
        disp(['Could not find variable that specifies cursor visibility.\n'
            'Make sure the variable is named correctly.'])
        keyboard
    end
else
    disp(['Could not find information about the start target in the TP table.\n',
        'Make sure the variable is named correctly.'])
end

%% Reach target position and size

% get column index
iCol = getColumnIndex(TPcolnames,'Reach_target'); % modify if target is named differently

% pre-allocate
nCol = length(iCol);
trialSpecs.xTarget = NaN(nTP,nCol);
trialSpecs.yTarget = NaN(nTP,nCol);
trialSpecs.targetRadius = NaN(nTP,nCol);
trialSpecs.targetWidth = NaN(nTP,nCol);
trialSpecs.targetHeight = NaN(nTP,nCol);

% get info from target table
for c = 1 : length(iCol)
    targetNumber = TPtable(:,iCol(c));
    trialSpecs.xTarget(notZero,c) = targetTable.X_GLOBAL(targetNumber);
    trialSpecs.yTarget(notZero,c) = targetTable.Y_GLOBAL(targetNumber);
    % radius
    if isfield(targetTable,'Radius')
        trialSpecs.targetRadius(notZero,c) = targetTable.Radius(targetNumber);
    elseif isfield(targetTable,'Visual_radius')
        trialSpecs.targetRadius(notZero,c) = targetTable.Visual_radius(targetNumber);
    end
    % horizontal and vertical size
    if isfield(targetTable,'Size_X')
        trialSpecs.targetWidth(notZero,c) = targetTable.Size_X(targetNumber);
        trialSpecs.targetHeight(notZero,c) = targetTable.Size_Y(targetNumber);
    end
end
% target shape
if any(~isnan(trialSpecs.targetWidth(:))) && any(trialSpecs.targetWidth(:)>0)
    trialSpecs.targetShape(1:nTP,1) = {'rect'};
elseif any(~isnan(trialSpecs.targetRadius(:))) && any(trialSpecs.targetRadius(:)>0)
    trialSpecs.targetShape(1:nTP,1) = {'circle'};
else
    disp('Check target size')
    keyboard
end
% haptic target position (identical to visual target position unless there
% is a perturbation to the cursor feedback)
trialSpecs.xTarget_haptic = trialSpecs.xTarget;
trialSpecs.yTarget_haptic = trialSpecs.yTarget;

%% Additional specs
% additional targets (from target table) could be added here in a similar way

if ~isempty(additionalSpecs)
    for i = 1 : length(additionalSpecs)
        % create field
        trialSpecs.(additionalSpecs{i}) = NaN(nTP,1);
        
        % get column index
        iCol = getColumnIndex(TPcolnames,additionalSpecs{i});
        
        % get info from TP table
        if iCol>0
            trialSpecs.(additionalSpecs{i})(notZero,1) = TPtable(:,iCol);
        end
    end
end
