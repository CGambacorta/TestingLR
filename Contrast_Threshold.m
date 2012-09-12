% Contrast_Threshold(fname, startContrast, noiseThresh, eccentricity)
% This program measures contrast threshold in the periphery using one
% staircase: 3 up/ 1 down (See staircase.m for more info.)
%
% IMPORTANT:
% To exit program while running:
% 1. Hit q during response time. If this doesn't work: 
% use the following key combinations:
% 1. ctrl-c                 (to allow MATLAB to recognize key inputs)
% 2. command/meta-period    (to stop the while loop)
% 3. sca                    (to close all screens)
%
% Inputs:
% fname = output file name.
% startContrast = starting contrast for the staircase.
% pset = parameter file, keep blank for default
%
% Example:
% Convolved_Threshold('TestData', 0.4, 0.1, 2.5)
% The start contrast level for the staircase will be 0.4, the stimulus
% eccentricity will be 2.5 degrees for the center and the result will be
% saved in 'TestData.mat'. (These are the default valuesif no argument is
% given.) Load TestData.mat and open the variable 'rec'to see all the
% information for every trial.
%
% Output:
% 'rec' is a 6 column variable which contains the following:
% 1st column: trial number
% 2nd column: location, 1(left),-1(right)
% 3rd column: up or down, 1(up), -1(down)
% 4th column: contrast of stimuli in this trial
% 5th column: reponse from subjects, 1(up), -1(down)
% 6th column: correct or not, 1(correct), 0 (wrong)
%
%
% Rachel Albert/ Christina Gambacorta (Levi Lab), 07/11/2012


function Contrast_Threshold(fname, startContrast, pset)

%------------------------SET PARAMETERS HERE---------------------------%
%File and Key Settings:
% This sets the default values if no arguments are given
if nargin<3 || isempty(pset), pset = 'default'; end
if nargin<2 || isempty(startContrast), startContrast = 0.4; end
if nargin<1 || isempty(fname), fname = 'TestData'; end

% Key codes
KbName('UnifyKeyNames');

% parameter set
[p w rect] = initialize(startContrast, pset);

ListenChar(2); % blocks input to matlab while function is running

%gray = 0.5; % TODO: check that this gray matches background color

%------------------------STAIRCASE AND OUTPUT--------------------------%
%Create Staircase
sc = staircase('create',[1 3],[0 1],[],p.nReversals);
sc.stimVal = startContrast;

%Order of Trials:
% left or right location, 1 or -1
loc = Shuffle([ones(p.stairTrials/2, 1); -ones(p.stairTrials/2, 1)]);

% 45 or -45 orientation, 1 or -1
ori = Shuffle([ones(p.stairTrials/2, 1); -ones(p.stairTrials/2, 1)]);

%Create 'rec' Variable:
% save trial number, location, and orientation in columns 1-3
rec = nan(p.stairTrials,7);
rec(:,1) = (1:p.stairTrials)';
rec(:,2) = loc;
rec(:,3) = ori;

%------------------------BEGIN MAIN FUNCTION---------------------------%
try
    % The Psychtoolbox AssertOpenGL command will issue an error message
    % if you try to execute this script on a computer without OpenGL.
    AssertOpenGL;

    p.fi = Screen('GetFlipInterval',w);   % fixation interval
    nf = round(p.stimDur/p.fi); % number of frames for gabor and noise
    p.stimDur = nf * p.fi; % update the real duration for record

    %Linearize Gamma: NOTE: This code may or may not work - use with caution!
    if p.linearize
        fid = fopen(p.gammaTableFile,'r');
        screen_clut = fread(fid,[256 3],'float64');
        fclose(fid);
        screen_clut = screen_clut/256;
        Screen('LoadNormalizedGammaTable',p.scrID,screen_clut);
    else
        Screen('ReadNormalizedGammaTable',p.scrID);
    end

    HideCursor;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Stimuli Setup:
    % fixation cross
    fix = ones(p.fixsz)*p.background/256;
    fix(p.fixsz/2 - p.fixwd/2 + 1: p.fixsz/2 + p.fixwd/2, :) = 0;
    fix(:, p.fixsz/2 - p.fixwd/2 + 1: p.fixsz/2 + p.fixwd/2) = 0;
    fixRect = CenterRect([0 0 p.fixsz p.fixsz],rect);
    fix=Screen('MakeTexture',w,fix,[],[],2);

    % generate gabor
    [x,y] = meshgrid(linspace(-p.radius, p.radius, p.halfPicSize * 2));
    gaborResult = exp(-((x.^2+y.^2)/2/p.sigma.^2)) .* ...
        sind(360*p.sf*(x*cosd(0) + y*sind(0))); % generates gabor

    % generate gaussian envelope
    gaussianSpaceConstant = p.sigma * p.ppd;
    [x, y] = meshgrid(-p.halfPicSize:p.halfPicSize - 1, -p.halfPicSize:p.halfPicSize - 1);
    envelope = exp(-((x .^ 2) + (y .^ 2)) /((sqrt(2) * gaussianSpaceConstant)^ 2));

    % generate noise pattern
    [fx, fy] = meshgrid(-p.halfPicSize:p.halfPicSize - 1, -p.halfPicSize:p.halfPicSize - 1);
    fx = fx./p.halfPicSize;
    fy = fy./p.halfPicSize;
    f = sqrt(fx.*fx+fy.*fy);
    f(p.halfPicSize+1, p.halfPicSize+1) = 1.0;
    lowcutoff = 0.05;
    highcutoff = 0.2;
    lowg = 1.0./(1.0+ (lowcutoff./f).^4);
    highg = 1.0./(1.0+(f./highcutoff).^4);
    g = lowg .* highg;
    g(p.halfPicSize+1,p.halfPicSize+1) = 1.0;
    Ntmp = sqrt(2) * erfinv(2 * rand(p.halfPicSize, p.halfPicSize) - 1);
    Ntmp = Expand(Ntmp, 2, 2);
    fNtmp = fftshift(fft2(Ntmp));
    fNtmp = fftshift(fNtmp.*g);
    Ntmp = real(ifft2(fNtmp));

    % add envelope to noise pattern
    noiseData = Ntmp;
    noiseResult = envelope.*noiseData;

    %Instructions 1:
    Screen(w,'TextSize',28); Screen(w,'TextFont','Arial');
    str1 = 'Please keep your eyes focused on the cross.\n\nUse your peripheral vision to judge the orientation of the pattern.\n\n';
    str2 = 'Press UP arrow for up-tiled and DOWN arrow for down-tilted.\n\n\n\n';
    str3 = 'Press any key to start.';
    str = strcat(str1,str2,str3);
    DrawFormattedText(w, str, 'center','center', 0);

    Screen(w,'Flip'); % show instruction
    KbPressWait(-3); % wait till any key to start
    Screen('Flip',w); % remove instruction


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Running the Experiment:
    p.start = datestr(now); % for record
    Priority(MaxPriority(w));

    i = 1; % trial index

    while ~isempty(sc.stimVal)
        if loc(i) == -1 %stim right
            stimrect = p.offsetRect(2, :);
            noiserect = p.offsetRect(1, :); %noise left
        elseif loc(i) == 1 %stim left
            noiserect = p.offsetRect(2, :); %noise right
            stimrect = p.offsetRect(1, :);
        end

        % create image matrixes
        %gabor = p.gray + gaborResult * sc.stimVal * p.gray;
        gabor = p.gray + (gaborResult * sc.stimVal * p.gray + noiseResult * p.noiseContrast * p.gray);
        %noiseMask = p.gray + noiseResult * p.noiseContrast * p.gray;
        noiseOnly = p.gray + (noiseResult * sc.stimVal * p.gray + noiseResult * p.noiseContrast * p.gray);
        %noisyGabor = gabor + noiseMask;
        noisetexture = Screen('MakeTexture',w,noiseOnly,[],[],2);
        stimtexture = Screen('MakeTexture',w,gabor,[],[],2);

        %Frames
            %slength = length(gaborResult);
            framexy = [-p.slength/2 p.slength/2 -p.slength/2 p.slength/2 ...
                -p.slength/2 -p.slength/2 p.slength/2 p.slength/2;...
                -p.slength/2 -p.slength/2 p.slength/2 p.slength/2 ...
                -p.slength/2 p.slength/2 -p.slength/2 p.slength/2];
        
%         %Frames
%         slength = length(gabor);
%         framexy = [-slength/2 slength/2 -slength/2 slength/2 ...
%             -slength/2 -slength/2 slength/2 slength/2;...
%             -slength/2 -slength/2 slength/2 slength/2 ...
%             -slength/2 slength/2 -slength/2 slength/2];
% 

% Present inital fixation cross
            Screen('DrawTexture',w,fix,[],p.fixRect,0,1); % display fix
            %Screen('DrawLines',w,framexy,1,0,p.srectc(1,:),0); %left frame
            %Screen('DrawLines',w,framexy,1,0,p.srectc(2,:),0); %right frame
            Screen('Flip',w);
        WaitSecs(p.soa  + .001*randi(400)); % add random wait

        % Present stimuli
        Screen('DrawTexture',w,fix,[],fixRect,0,1);
        Screen('DrawTexture',w,stimtexture,[],stimrect,p.dAngle*ori(i),1);
        Screen('DrawTexture',w,noisetexture,[],noiserect,p.dAngle*ori(i),1);
        %Screen('DrawLines',w,framexy,1,0,p.srectc(1,:),0); %left frame
        %Screen('DrawLines',w,framexy,1,0,p.srectc(2,:),0); %right frame
        Screen('Flip',w);

        WaitSecs(p.stimDur);

        % display fixation
        Screen('DrawTexture',w,fix,[],p.fixRect,0,1);
        %Screen('DrawLines',w,framexy,1,0,p.srectc(1,:),0); %left frame
        %Screen('DrawLines',w,framexy,1,0,p.srectc(2,:),0); %right frame
        Screen('Flip',w);

        % get response
       tic
       response = 0; 
       
       while toc < p.respDur;
           [down time kcode] = KbCheck(-3);
           if down
               KeyPress = KbName(kcode);
               if strcmp(KeyPress,'UpArrow')
                   response = 1;
                   rec(i,5) = response;
                   %rec(i,6) = (response == rec(i,5));
                   rec(i,7) = toc;
               elseif strcmp(KeyPress,'DownArrow')
                   response = -1;
                   rec(i,5) = response;
                   %rec(i,6) = (response == rec(i,5));
                   rec(i,7) = toc;
               elseif strcmp(KeyPress,'q')
                   error('esc:exit','ESC pressed.');
               elseif strcmp(KeyPress, ~'DownArrow') || strcmp(KeyPress, ~'UpArrow') || strcmp(KeyPress, ~'q')
                   down = 0;
               end
           end
       end

        % At the end of each trial save variables in rec
        %rec(i,5) = response;
        
        if abs(rec(i,5)) < 1 % if no response, not correct
            data = 0;
        else
        data = response == rec(i,3); % compare response
        end
        
        rec(i,4) = sc.stimVal; % record contrast value
        rec(i,6) = data;
        
        cd Contrast_Threshold_Data % save in data folder
        save(fname,'p','sc','rec');
        cd ..
        sc = staircase('update',sc,data); % update staircase
        clear stim;
        i = i + 1; % update trial number

        % Feedback sounds
        if data == 1
            sound(p.beep);
        elseif data == 0
            sound(p.beepError);
        end
    end

    cd Contrast_Threshold_Data % save in data folder
    save(fname,'p','sc','rec');
    cd ..

    % Compute threhold with discarded reversals
    sc = staircase('compute',sc, p.nDiscard);

    % Plot and save staircase figure
    staircase('plot', sc);
    dataPlot = figure(7);
    cd Contrast_Threshold_Data
    ts = DATESTR(now,1);
    filename = sprintf('%s_staircase_%s_.fig', fname, ts);
    saveas(dataPlot, filename ,'fig')
    cd ..

    
catch whatiserror % in case of error or user exit
    ShowCursor;
    Priority(0);
    ListenChar(0);
    Screen('Closeall');
    whatiserror.message
    whatiserror.stack
end

Screen('CloseAll');
Priority(0);
ListenChar(0);

end % main function ends

