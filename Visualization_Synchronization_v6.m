%% =========================================================
% EEG + EOG + TANGENTIAL IMU + VIDEO 
% Developed by: Lianne Sanchez-Rodriguez & Maxine Annel Pacheco-Ramirez
% Affiliation: University of Houston IUCRC BRAIN Center
% Date: September 4, 2026
%% =========================================================
clear; clc; close all;

% R2019b-compatible stability setting
opengl software

%% FILES
eegFile   = 'C:\Users\lsanch46\Downloads\Chicago_Session_split\participant_01.set';
imuFile   = 'C:\Users\lsanch46\Downloads\Opals_Post_Chicago_Split\participant_01_post_sensor_742.mat';
videoFile = 'C:\Users\lsanch46\Downloads\Chicago_Video.mp4';

%% SETTINGS
FsTarget = 100;
FsIMU    = 128;
t0Video  = 0;
winSec   = 5;

% move these channels to EOG panel
eogCh = [15 20 31 32];

tangentialAxis   = 1;   % 1=X,2=Y,3=Z
useNormalizedIMU = false;
imuYLim          = [-1 1];   % increased for better visibility
eegYLim          = [];
eogYLim          = [];
eegSpacing       = 4;
eegLineWidth     = 0.6;

trigColor        = [1 0 1];
trigLineWidth    = 2.0;
showTriggerText  = true;

% performance controls
plotUpdateHz   = 15;    % EEG/EOG/IMU redraw rate
videoSkip      = 1;     % 1=no skip, 2=every other frame
eegDecimPlot   = 2;     % plot fewer points each update

%% EEGLAB
if exist('pop_loadset','file') ~= 2
    error('pop_loadset not found. Run eeglab once.');
end
EEG = pop_loadset(eegFile);

% Close EEGLAB GUI if open (reduce graphics load)
f = findall(0,'Type','figure','Name','EEGLAB');
if ~isempty(f), close(f); end

%% EEG -> 100 Hz
FsEEG_orig = double(EEG.srate);
eegRawOrig = double(EEG.data); % nChan x nSamples
[nChan,~]  = size(eegRawOrig);

[pE,qE] = rat(FsTarget/FsEEG_orig,1e-12);
eegRaw = resample(eegRawOrig',pE,qE)'; % nChan x nEEG
FsEEG = FsTarget;
[nChan,nEEG] = size(eegRaw);
tEEG = (0:nEEG-1)/FsEEG;

%% Split EEG main vs EOG channels
allCh = 1:nChan;
eogCh = eogCh(eogCh>=1 & eogCh<=nChan);
eegKeepCh = setdiff(allCh, eogCh, 'stable');

eegMainRaw = eegRaw(eegKeepCh,:); % channels for EEG panel
eogRaw     = eegRaw(eogCh,:);     % channels for EOG panel

nMain = numel(eegKeepCh);
nEOG  = numel(eogCh);

%% IMU load -> 100 Hz
M = load(imuFile);
fns = fieldnames(M);
imu = [];
for i = 1:numel(fns)
    v = M.(fns{i});
    if isnumeric(v) && ismatrix(v) && (size(v,2)==3 || size(v,1)==3)
        imu = double(v);
        fprintf('Using IMU variable: %s\n', fns{i});
        break;
    end
end
if isempty(imu), error('No Nx3/3xN IMU matrix found.'); end
if size(imu,2)~=3, imu = imu'; end

[pI,qI] = rat(FsTarget/FsIMU,1e-12);
imuResCell = cell(1,3);
nResEach = zeros(1,3);
for k = 1:3
    imuResCell{k} = resample(imu(:,k),pI,qI);
    nResEach(k) = numel(imuResCell{k});
end
nRes = min(nResEach);
imu100 = zeros(nRes,3);
for k = 1:3
    imu100(:,k) = imuResCell{k}(1:nRes);
end

% align IMU length to EEG length
if nRes >= nEEG
    imuOnEEG = imu100(1:nEEG,:);
else
    imuOnEEG = [imu100; repmat(imu100(end,:), nEEG-nRes, 1)];
end
imuTan = imuOnEEG(:,tangentialAxis);

%% VIDEO + OVERLAP
if ~isfile(videoFile), error('Video not found: %s', videoFile); end
vrInfo = VideoReader(videoFile);

tEnd = min([tEEG(end), t0Video + vrInfo.Duration]);
keep = tEEG <= tEnd;

t = tEEG(keep);
eegMainRaw = eegMainRaw(:,keep);
eogRaw     = eogRaw(:,keep);
imuTan     = imuTan(keep);

%% PREP EEG main
eegZ = eegMainRaw - mean(eegMainRaw,2);
sd = std(eegZ,0,2); sd(sd==0)=1;
eegZ = eegZ./sd;
offsets = ((nMain-1):-1:0)'*eegSpacing;

labels = "Ch " + string(eegKeepCh);
if isfield(EEG,'chanlocs') && ~isempty(EEG.chanlocs)
    for i = 1:nMain
        c = eegKeepCh(i);
        if c <= numel(EEG.chanlocs) && isfield(EEG.chanlocs(c),'labels') && ~isempty(EEG.chanlocs(c).labels)
            labels(i)=string(EEG.chanlocs(c).labels);
        end
    end
end
if isempty(eegYLim)
    eegYLim = [min(offsets)-2.5, max(offsets)+2.5];
end

%% PREP EOG (same scaling style as EEG for comparability)
eogZ = eogRaw - mean(eogRaw,2);
sdE = std(eogZ,0,2); sdE(sdE==0)=1;
eogZ = eogZ./sdE;

eogSpacing = eegSpacing; % same spacing as EEG
eogOffsets = ((nEOG-1):-1:0)'*eogSpacing;

eogLabels = "Ch " + string(eogCh);
if isfield(EEG,'chanlocs') && ~isempty(EEG.chanlocs)
    for i = 1:nEOG
        c = eogCh(i);
        if c <= numel(EEG.chanlocs) && isfield(EEG.chanlocs(c),'labels') && ~isempty(EEG.chanlocs(c).labels)
            eogLabels(i)=string(EEG.chanlocs(c).labels);
        end
    end
end
if isempty(eogYLim)
    eogYLim = [min(eogOffsets)-2.5, max(eogOffsets)+2.5];
end

%% PREP IMU
if useNormalizedIMU
    imuTanPlot = imuTan - mean(imuTan);
    s = std(imuTanPlot); if s==0, s=1; end
    imuTanPlot = imuTanPlot/s;
else
    imuTanPlot = imuTan;
end

%% TRIGGERS
evtT = [];
evtLbl = strings(0,1);
if isfield(EEG,'event') && ~isempty(EEG.event)
    for k = 1:numel(EEG.event)
        if isfield(EEG.event(k),'latency') && ~isempty(EEG.event(k).latency)
            tk = (double(EEG.event(k).latency)-1)/FsEEG_orig;
            if tk>=t(1) && tk<=t(end)
                evtT(end+1,1)=tk; %#ok<AGROW>
                if isfield(EEG.event(k),'type') && ~isempty(EEG.event(k).type)
                    evtLbl(end+1,1)=string(EEG.event(k).type); %#ok<AGROW>
                else
                    evtLbl(end+1,1)="evt"; %#ok<AGROW>
                end
            end
        end
    end
end

%% FIGURE
fig = figure('Color','w','Units','normalized','OuterPosition',[0.02 0.05 0.96 0.9], ...
    'Name','SPACE pause/resume | Q/ESC quit');
set(fig,'CurrentCharacter',char(0),'Renderer','opengl');

% Left video, right stacked panels: EEG top, EOG middle, IMU thin bottom
axV = axes('Parent',fig,'Position',[0.03 0.08 0.46 0.84]); % video = 50% width
ax1 = axes('Parent',fig,'Position',[0.56 0.48 0.41 0.44]); % EEG taller
ax3 = axes('Parent',fig,'Position',[0.56 0.22 0.41 0.20]); % EOG smaller
ax2 = axes('Parent',fig,'Position',[0.56 0.08 0.41 0.08]); % IMU thin

% video image created once
vr = VideoReader(videoFile);
if hasFrame(vr), fr = readFrame(vr); else, fr = zeros(480,640,3,'uint8'); end
imH = image(axV,fr);
axis(axV,'image'); axis(axV,'off'); title(axV,'Video');

% EEG panel
hold(ax1,'on');
hEEG = gobjects(nMain,1);
for c = 1:nMain
    hEEG(c) = plot(ax1,nan,nan,'k','LineWidth',eegLineWidth);
end
yticks(ax1,sort(offsets));
yticklabels(ax1,labels(end:-1:1));
ylabel(ax1,'EEG');
%title(ax1,'EEG stacked');
grid(ax1,'on');
cursor1 = xline(ax1,0,'--r','LineWidth',1.2);
ylim(ax1,eegYLim);
ax1.XTick = [];

% EOG panel
hold(ax3,'on');
hEOG = gobjects(nEOG,1);
for c = 1:nEOG
    hEOG(c) = plot(ax3,nan,nan,'b','LineWidth',1.0);
end
yticks(ax3,sort(eogOffsets));
yticklabels(ax3,eogLabels(end:-1:1));
ylabel(ax3,'EOG');
%title(ax3,'EOG');
grid(ax3,'on');
cursor3 = xline(ax3,0,'--r','LineWidth',1.2);
ylim(ax3,eogYLim);
ax3.XTick = [];

% IMU panel
hold(ax2,'on');
hTan = plot(ax2,nan,nan,'r','LineWidth',1.4);
ylabel(ax2,'IMU'); xlabel(ax2,'Time (s)');
%title(ax2,'IMU Tangential');
grid(ax2,'on');
cursor2 = xline(ax2,0,'--r','LineWidth',1.2);
ylim(ax2,imuYLim);

% Link x-axes
linkaxes([ax1 ax3 ax2],'x');
xlim(ax1,[0 min(winSec,t(end))]);

% Trigger lines (keep on all panels)
for k = 1:numel(evtT)
    xline(ax1,evtT(k),':','Color',trigColor,'LineWidth',trigLineWidth,'HandleVisibility','off');
    xline(ax3,evtT(k),':','Color',trigColor,'LineWidth',trigLineWidth,'HandleVisibility','off');
    xline(ax2,evtT(k),':','Color',trigColor,'LineWidth',trigLineWidth,'HandleVisibility','off');
end

% Trigger label ONLY on EEG panel
if showTriggerText
    hTrigTxt1 = text(ax1,0.00,1.04,'','Units','normalized', ...
        'VerticalAlignment','bottom','HorizontalAlignment','left', ...
        'Color',trigColor,'FontWeight','bold','BackgroundColor','w','Margin',2,'Clipping','off');
else
    hTrigTxt1 = gobjects(1);
end

%% PLAYBACK
isPaused = false; stopNow = false;
fps = vr.FrameRate; dtFrame = 1/fps;
nFrames = max(1,floor(vr.Duration*fps));
frameIdx = 2;  % first frame already displayed

tWall0 = tic; pauseAccum = 0; pauseTic = [];
lastPlotUpdate = -inf; plotDt = 1/plotUpdateHz;

while frameIdx <= nFrames
    if ~ishandle(fig), break; end

    % keyboard controls
 ch = get(fig,'CurrentCharacter');
if ~isempty(ch) && ch~=char(0)
    % Quit on Q, ESC, ENTER, or SPACE
    if ch=='q' || ch=='Q' || double(ch)==27 || ch==char(13) || ch==char(10) || ch==' '
        stopNow = true;
    end
    set(fig,'CurrentCharacter',char(0));
end
    if stopNow, break; end

    if isPaused
        drawnow;
        pause(0.01);
        continue;
    end

    tExpected = (frameIdx-1)*dtFrame;
    tElapsed = toc(tWall0)-pauseAccum;

    if tElapsed < tExpected
        pause(min(0.003,tExpected-tElapsed));
        drawnow;
        continue;
    end

    % video update
    if mod(frameIdx-1,videoSkip)==0 && hasFrame(vr)
        fr = readFrame(vr);
        set(imH,'CData',fr);
    elseif hasFrame(vr)
        readFrame(vr); % keep timing
    end

    tNow = tExpected + t0Video;
    if tNow > t(end), break; end

    % throttle signal redraw
    if (tNow-lastPlotUpdate) >= plotDt
        t0 = max(0,tNow-winSec);
        t1 = min(t0+winSec,t(end));
        idx = find(t>=t0 & t<=t1);

        if ~isempty(idx)
            idx = idx(1:eegDecimPlot:end);
            tw = t(idx);

            % EEG main
            for c = 1:nMain
                set(hEEG(c),'XData',tw,'YData',eegZ(c,idx)+offsets(c));
            end

            % EOG
            for c = 1:nEOG
                set(hEOG(c),'XData',tw,'YData',eogZ(c,idx)+eogOffsets(c));
            end

            % IMU
            set(hTan,'XData',tw,'YData',imuTanPlot(idx));

            xlim(ax1,[t0 t1]); xlim(ax3,[t0 t1]); xlim(ax2,[t0 t1]);
            cursor1.Value = tNow; cursor3.Value = tNow; cursor2.Value = tNow;

            % Trigger text ONLY on EEG panel
            if showTriggerText
                if ~isempty(evtT)
                    iPast = find(evtT<=tNow,1,'last');
                    if ~isempty(iPast)
                        msg = sprintf('Trigger: %s @ %.3fs', char(evtLbl(iPast)), evtT(iPast));
                    else
                        msg = 'Trigger: none yet';
                    end
                else
                    msg = 'Trigger: none';
                end
                set(hTrigTxt1,'String',msg);
            end
        end
        lastPlotUpdate = tNow;
    end

    drawnow limitrate nocallbacks;
    %drawnow nocallbacks;
    frameIdx = frameIdx + 1;
end

disp('Done.');