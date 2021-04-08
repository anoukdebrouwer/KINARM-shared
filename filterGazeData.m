function [xyGaze,vxyGaze,blink,notUpdated] = filterGazeData(xyGaze_raw,xyLimits,fs,nsMax,tGaze)
% filterGazeData  Remove blinks and filters raw gaze data from Eyelink eye tracker
%
% xyGaze = filterGazeData(xyGaze_raw,xyLimits,fs) filters the x and y gaze
% positions in xyGaze_raw obtained from an Eyelink eye tracker at sampling
% frequency fs (typically 500 Hz) and returns them in output xyGaze. Blinks
% or intervals in which the pupil was lost (i.e., outside the area
% [xmin xmax ymin ymax] specified in xyLimits) for more than 20 samples are
% removed, shorter intervals are linearly interpolated. The x and y gaze
% positions are low-pass filtered with a 2nd order Butterworth filter with
% 30 Hz cutoff frequency. Arrays with NaNs are returned if the trial did
% not contain valid gaze data.
%
% xyGaze = filterGazeData(xyGaze_raw,xyLimits,fs,nsMax) interpolates data
% losses up to nsMax samples, and classifies intervals longer than nsMax
% samples as blinks.
%
% xyGaze = filterGazeData(xyGaze_raw,xyLimits,fs,nsMax,tGaze) resamples the
% data to regular sampling if sampling is irregular, based on the time
% stamps in tGaze.
%
% [xyGaze,vxyGaze] = filterGazeData(__) returns the x and y gaze velocity
% in output vxyGaze corresponding to the positions in xyGaze.
%
% [xyGaze,vxyGaze,blink] = filterGazeData(__) additionally returns a
% boolean variable in output blink that is true during blinks and false
% otherwise.
%
% [xyGaze,vxyGaze,blink,notUpdated] = filterGazeData(__) additionally
% returns a boolean variable in output notUpdated that indicates where
% sampling problems occured (no new gaze position for at least 20 samples).

% MIT License
% Copyright (c) 2020 Anouk de Brouwer

plotFilterResults = false; % set to true to plot raw and filtered data

if isstruct(xyLimits) % compatible with older code
    Exp = xyLimits;
    fs = Exp.setup.fs;
    xyLimits = [-1.5 1.5 -1.5 1.5]*Exp.stim.targetDistance;
    nsMax = 20;
end
maxVel = 1500; % cm/s or deg/s

nSamples = length(xyGaze_raw);

% set defaults if not provided
resampleGaze = false;
if nargin==3
    nsMax = 20;
    tGaze = [0:1/fs:(nSamples-1)/fs]';
elseif nargin==2 || nargin==4
    tGaze = [0:1/fs:(nSamples-1)/fs]';
else
    tGaze = tGaze - tGaze(1);
    resampleGaze = true;
end

% compute gaze velocity from raw data - for plot only
vxyGaze_raw(:,1) = gradient(xyGaze_raw(:,1),1/fs);
vxyGaze_raw(:,2) = gradient(xyGaze_raw(:,2),1/fs);

%% Resample gaze position if necessary

% resample if original sampling is irregular
if resampleGaze
    t = [0:1/fs:(nSamples-1)/fs]';
    nw = [true; tGaze(2:end)~=tGaze(1:end-1)]; % new sample
    x = interp1(tGaze(nw),xyGaze_raw(nw,1),t);
    y = interp1(tGaze(nw),xyGaze_raw(nw,2),t);
    if isnan(x(end))
        xnan = isnan(x);
        x(xnan) = x(find(~xnan,1,'last'));
        y(xnan) = y(find(~xnan,1,'last'));
    end
    % for plot only
    vx_raw = gradient(xyGaze_raw(nw,1),1/fs);
    vy_raw = gradient(xyGaze_raw(nw,2),1/fs);
    vx = interp1(tGaze(nw),vx_raw,t);
    vy = interp1(tGaze(nw),vy_raw,t);
    vxyGaze_raw = [vx vy];
    % update
    xyGaze_raw = [x y];
end

%% Detect blinks and interpolate sampling issues

xyGaze_int = xyGaze_raw;
blink = false(nSamples,1);
interpGaze = false(nSamples,1);

% detect blinks and sampling issues [(x,y) is (-100,-100)]
% intervals in which gaze is outside workspace
% or in which there is a large difference in gaze position
dxyGaze_raw = diff(xyGaze_raw);
dGaze_raw = [0; sqrt(dxyGaze_raw(:,1).^2 + dxyGaze_raw(:,2).^2)];
gazeOut = xyGaze_raw(:,1)<xyLimits(1) | xyGaze_raw(:,1)>xyLimits(2) | ...
    xyGaze_raw(:,2)<xyLimits(3) | xyGaze_raw(:,2)>xyLimits(4) | ...
    dGaze_raw>(maxVel/fs);
[~,iStartEnd] = findIntervals(gazeOut,1);

% if there are blinks and/or sampling issues, get duration and precompute padding
if ~isnan(iStartEnd)
    % combine intervals with issues if there are less than 4 samples (= 2 new
    % samples at 500 Hz) between them
    if size(iStartEnd,1)>1
        d = iStartEnd(2:end,1)-iStartEnd(1:end-1,2);
        shortd = find(d<4,1);
        while ~isempty(shortd)
            iStartEnd(shortd+1,1) = iStartEnd(shortd,1);
            iStartEnd(shortd,:) = [];
            if size(iStartEnd,1)>1
                d = iStartEnd(2:end,1)-iStartEnd(1:end-1,2);
                shortd = find(d<4,1);
            else
                shortd = [];
            end
        end
    end
    % compute padding interval
    % should blink padding be different for tablet and kinarm?
    n = iStartEnd(:,2)-iStartEnd(:,1)+1;
    %iPaddedBlink = [iStartEnd(:,1)-10 iStartEnd(:,2)+44]; % tablet
    iPaddedBlink = [iStartEnd(:,1)-30 iStartEnd(:,2)+100]; % kinarm
    iPaddedBlink(iPaddedBlink<1) = 1;
    iPaddedBlink(iPaddedBlink>nSamples) = nSamples;
    
    % loop over blinks/sampling issues
    for i = 1 : size(iStartEnd,1)
        % interpolate sampling issues where number of samples <= nsMax
        if n(i) <= nsMax
            iVal = [iStartEnd(i,1)-1 iStartEnd(i,2)+1];
            % do not interpolate if signal loss occurs at start or end of recording
            if iVal(1)<1 || iVal(2)>nSamples
                blink(iStartEnd(i,1):iStartEnd(i,2)) = true;
            else % interpolate
                x_int = interp1(iVal,xyGaze_int(iVal,1),iVal(1):iVal(end),'linear');
                y_int = interp1(iVal,xyGaze_int(iVal,2),iVal(1):iVal(end),'linear');
                % only interpolate within limits
                while any(x_int<xyLimits(1)) || any(x_int>xyLimits(2)) || ...
                        any(y_int<xyLimits(3)) || any(y_int>xyLimits(4))
                    iVal(end) = iVal(end)+1;
                    x_int = interp1(iVal,xyGaze_int(iVal,1),iVal(1):iVal(end),'linear');
                    y_int = interp1(iVal,xyGaze_int(iVal,2),iVal(1):iVal(end),'linear');
                end
                xyGaze_int(iVal(1):iVal(end),:) = [x_int' y_int'];
                interpGaze(iVal(1):iVal(end),1) = true;
            end
        else % pad if real blink (n>nsMax)
            blink(iPaddedBlink(i,1):iPaddedBlink(i,2)) = true;
        end
    end
end

%% Check sampling of gaze position

notUpdated = false(nSamples,1);
dx = [0; diff(xyGaze_raw(:,1))];
dy = [0; diff(xyGaze_raw(:,2))];
[~,iNotUpdated] = findIntervals([dx==0 dy==0 ~gazeOut],20); % should be updated every sample
if ~isnan(iNotUpdated)
    % save non-updated samples
    notUpdated(iNotUpdated(1,1):iNotUpdated(1,2)) = true;
end
blink(notUpdated) = false; % make sure the padding doesn't interfere with not updated samples

%% Apply low-pass butterworth filter and compute gaze velocity

noData = blink | notUpdated;
[~,iYesData] = findIntervals(~noData,7); % at least 6 samples needed to filter the data

if ~isnan(iYesData)
    
    [b,a] = butter(2,30/(fs/2));
    xyGaze = NaN(size(xyGaze_int));
    vxyGaze = NaN(size(xyGaze_int));
    for i = 1 : size(iYesData,1);
        int = iYesData(i,1):iYesData(i,2);
        % filter gaze position (excluding blinks)
        xyGaze(int,:) = filtfilt(b,a,xyGaze_int(int,:));
        % compute velocity (excluding blinks)
        vxyGaze(int,1) = gradient(xyGaze(int,1),1/fs);
        vxyGaze(int,2) = gradient(xyGaze(int,2),1/fs);
    end
    
    %% Plot to check
    
    % check removing of blinks and filtering
    if plotFilterResults && ~isnan(iStartEnd(1))
        % x position
        subplot(2,2,1)
        plot([xyGaze_raw(:,1),xyGaze(:,1)],'.-')
        legend('raw','blinks removed and filtered')
        axis([1 nSamples xyLimits(1:2)])
        xlabel('Time (samples)')
        ylabel('x')
        % y position
        subplot(2,2,2)
        plot([xyGaze_raw(:,2),xyGaze(:,2)],'.-')
        axis([1 nSamples xyLimits(3:4)])
        xlabel('Time (samples)')
        ylabel('y')
        % x velocity
        subplot(2,2,3)
        plot([vxyGaze_raw(:,1),vxyGaze(:,1)],'.-')
        legend('raw','from filtered position data')
        axis([1 nSamples -1000 1000])
        xlabel('Time (samples)')
        ylabel('vx')
        % y velocity
        subplot(2,2,4)
        plot([vxyGaze_raw(:,2),vxyGaze(:,2)],'.-')
        axis([1 nSamples -1000 1000])
        xlabel('Time (samples)')
        ylabel('vy')
        keyboard
        %waitforbuttonpress
    end
    
else
    
    %% No gaze data in trial
    
    xyGaze = NaN(size(xyGaze_raw));
    vxyGaze = NaN(size(xyGaze_raw));
    
end
