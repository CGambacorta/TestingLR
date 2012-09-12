function [p w rect] = initialize(contrast, pset)
if strcmp(pset, 'default') ~= 1
    [p w rect] = optionalParameters(contrast, pset);
else
%------------------------SET PARAMETERS HERE---------------------------%
%Experimental Paremters:
% monitor related parameters
p.distance          = 60;                        % cm, screen-eye distance
p.screenHeight      = 30.5;                      % cm, height of display region
p.scrID             = 0; %max(Screen('Screens'));    % which screen to display, set to 0 if running on laptop
p.linearize         = 0;                         % 0 = non-linearize; 1 = linearize
p.gammaTableFile    = 'Gamma_Table';             % Gamma Table File

% time parameters
p.soa               = 0.8;                       % duration of fixation cross
p.soaJitter         = 400;                       % max ms of jitter
p.iti               = 1;                         % interval from resp to next trial
p.waitDur           = 0.5;                       % interval from stim to response
p.stimDur           = 0.05;    %0.1;             % gabor display time
p.cueDur            = 0.1;                       % cue display time
p.respDur            = 1;                         % max response time

% stimulus parameters
p.cueType           = 1;                         % cue type; 0 = none; 1 = audio; 2 = visual
p.background        = 128;                       % grey value for background
p.fixwd             = 4;                         % fixation thickness in pixels
p.fixsz             = 24;                        % pixels, fixation cross size
p.radius            = 0.75;  %1.25;              % radius of stimulus in degrees
p.eccentricity      = 5;      %8;                % degree from fix to gabor center
p.sf               = 3;     %3.5;                % cycles/degree of gabor
p.sigma             = 0.25;    %0.5;             % degree, gaussian sigma (envelope)
p.lowCutoff         = 0.05;                      % noise patch low cutoff freq
p.highCutoff        = 0.2;                       % noise patch high cutoff freq
p.dAngle            = -95;                       % orientation of gabor (-95 is up, 95 is down)
p.gaborContrast     = contrast *1.15;           % gabor contrast
p.noiseContrast     = 0.33;                      % noise contrast
p.nBlocks           = 3;                         % Number of blocks
p.trialsPerBlock    = 40;                        % Number of trials per block
p.nTrials           = p.nBlocks * p.trialsPerBlock;  % number of trials

%Screen Setup:
Screen('Preference', 'VisualDebugLevel', 0);
PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'AllViews', 'EnableCLUTMapping');
[w rect] = Screen('OpenWindow', p.scrID, p.background);

% exact timing parameters
p.fixInt = Screen('GetFlipInterval', w);        % fixation interval
nFrames = round(p.stimDur/p.fixInt);            % #frames for stim
p.stimDur = nFrames * p.fixInt;                 % update for record

% exact size parameters
p.ppd = p.distance/p.screenHeight * tand(1) * rect(4);  % pixels per deg
p.pps = 2 * p.radius * p.ppd;                           % pixels per stim
p.halfPicSize = round(p.radius * p.ppd);                % 1/2 length of stim

% fix location
p.fixRect = CenterRect([0 0 p.fixsz p.fixsz], rect);

% stimulus location
p.crect = CenterRect([0 0 p.pps p.pps], rect);
p.offset = round(p.eccentricity/sqrt(2) * p.ppd); % stim offset from center
p.offsetRect = [OffsetRect(p.crect, -p.offset, p.offset); OffsetRect(p.crect, p.offset, p.offset)];

% for lateral locations
%p.offset = round(p.eccentricity * p.ppd); % stim offset from center lateral
%p.offsetRect = [OffsetRect(p.crect,-p.offset,0); OffsetRect(p.crect,p.offset,0)]; % lateral

% generate gabor
[x,y] = meshgrid(linspace(-p.radius, p.radius, p.halfPicSize * 2));
gaborResult = exp(-((x.^2+y.^2)/2/p.sigma.^2)) .* ...
                sind(360*p.sf*(x*cosd(0) + y*sind(0))); % generates gabor matrix

% frames locations
p.slength = length(gaborResult);
p.srect = CenterRect([0 0 p.slength p.slength], rect);     
p.srect = [OffsetRect(p.srect, -p.offset, p.offset);OffsetRect(p.srect, p.offset, p.offset)];
p.srectc = [(p.srect(1,3)-p.srect(1,1))/2+p.srect(1,1), (p.srect(1,4)-p.srect(1,2))/2+p.srect(1,2);...
                (p.srect(2,3)-p.srect(2,1))/2+p.srect(2,1),(p.srect(2,4)-p.srect(2,2))/2+p.srect(2,2)];
         
% for lateral locations            
%p.srect = [OffsetRect(p.srect, -p.offset, 0);OffsetRect(p.srect, p.offset, 0)];
%p.srectc = [(p.srect(1,3)-p.srect(1,1))/2+p.srect(1,1), (p.srect(1,4)-p.srect(1,2))/2+p.srect(1,2);...
                %(p.srect(2,3)-p.srect(2,1))/2+p.srect(2,1), (p.srect(2,4)-p.srect(2,2))/2+p.srect(2,2)];

% exact color definitions
% black = BlackIndex(p.scrID);
% white = WhiteIndex(p.scrID);
% TODO: Use this code?
p.gray = 0.5;

% sound parameters
p.beep      = sin(1:0.5:100); % correct sound
p.beepError = sin(1:0.25:100); % error sound
p.leftCue   = wavread('left_short_loud.wav');
p.rightCue  = wavread('right_short_loud.wav');

% staircase parameters
p.stairTrials =  200;                     % total trials
p.nReversals = 14;                    % number of reversals
p.nDiscard = 7;                       % num. initial staircases discarded
end
end

function [p w rect] = optionalParameters(contrast, pset)
[p w rect] = initialize(contrast, 'default');
if strcmp(pset, 'test') == 1
    p.stimDur = 1;
    p.gaborContrast = 0.5;
    p.nBlocks = 2;
    p.trialsPerBlock = 32;
    p.nTrials = p.nBlocks * p.trialsPerBlock;
end
end
    







