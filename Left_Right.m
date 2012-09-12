% Left_Right(fname, contrast, noiseThresh, eccentricity)
% This program measures accuracy and response time for a gabor
% tilt-discrimination task with valid and invalid audio cues
%
% IMPORTANT:
% To exit program while running:
% 1. Hit q during response time. If that doesn't work try the following:
% 1. ctrl-c                 (to allow MATLAB to recognize key inputs)
% 2. command/meta-period    (to stop the while loop)
% 3. sca                    (to close all screens)
%
% Inputs:
% fname = output file name.
% startContrast = contrast for the gabor stimulus
% pset = parameters, use blank for default.
%
% Example:
% CTProgram('subject', 0.2, 0.1, 2.5)
% The gabor contrast will be 0.2, the noise contrast will 0.1 * p.n,
% and the stimuli will be scaled for 2.5 degrees of eccentricity.
% The result will be saved in 'subject.mat'. (These are the default values
% if no argument is given.) Load subject.mat and open the variable 'rec'
% to see all the information for every trial.
%
% Output:
% 'rec' is a 6 column variable which contains the following:
% 1st column: trial number
% 2nd column: trial type (total of 10 types)
% 3rd column: audio cue, 1(left),-1(right)
% 4th column: location, 1(left),-1(right)
% 5th column: orientation up or down, 1(up), -1(down)
% 6th column: reponse from subjects, 1(up), -1(down)
% 7th column: correct or not, 1(correct), 0 (wrong)
% 8th column: response time
%
%
% Rachel Albert/ Christina Gambacorta (Levi Lab), 07/10/2012
% EDITED 8/24/12 to add % correct after each block, new proportions, no noise response.
% Also changed stimulus to be additive (not multiplicative).

function Left_Right(fname,contrast, pset)

%------------------------SET PARAMETERS HERE---------------------------%
%File and Key Settings:
% This sets the default values if no arguments are given
if nargin<3 || isempty(pset), pset = 'default'; end
if nargin<2 || isempty(contrast), contrast = 0.1; end
if nargin<1 || isempty(fname), fname = 'Subject'; end

% parameter set
[p w rect] = initialize(contrast, pset);
p.fname = fname;

% Key codes
KbName('UnifyKeyNames');

%---------------------TRIAL ORDER AND DEFINITIONS-----------------------%

% Headers for trialType: Type, Cue, Location, Orientation
% 40 trials per block (or can double to have 80), this is what determines
% the ratio. 

trialType = ...
    [1, -1, 1,  1; % Invalid, left, up 
    1, -1,  1,  1;
    2, -1,  1, -1; % Invalid, left, down 
    2, -1,  1, -1;
    3,  1, -1,  1; % Invalid, right, up 
    3,  1, -1,  1;
    4,  1, -1, -1; % Invalid, right, down 
    4,  1, -1, -1;
    5,  1,  1,  1; % Valid, left, up 
    5,  1,  1,  1;
    5,  1,  1,  1;
    5,  1,  1,  1;
    5,  1,  1,  1;
    5,  1,  1,  1;
    6,  1,  1, -1; % Valid, left, down 
    6,  1,  1, -1;
    6,  1,  1, -1;
    6,  1,  1, -1;
    6,  1,  1, -1;
    6,  1,  1, -1;
    7, -1, -1,  1; % Valid, right, down 
    7, -1, -1,  1;
    7, -1, -1,  1;
    7, -1, -1,  1;
    7, -1, -1,  1;
    7, -1, -1,  1;
    8, -1, -1, -1; % Valid, right, down 
    8, -1, -1, -1;
    8, -1, -1, -1;
    8, -1, -1, -1;
    8, -1, -1, -1;
    8, -1, -1, -1;
    9,  1,  0, -1; % Left cue, Noise, noise 
    9,  1,  0,  1;
    9,  1,  0, -1;
    9,  1,  0, -1;
    10, -1, 0, -1; % Right cue, Noise, noise 
    10, -1, 0, -1;
    10, -1, 0, -1;
    10, -1, 0,  1];

trialTypeSize = size(trialType,1);
assert((mod(p.nTrials,trialTypeSize) == 0), sprintf('Number of Trials must be multiple of %d',trialTypeSize));


%------------------------BEGIN MAIN FUNCTION---------------------------%
try
    % The Psychtoolbox AssertOpenGL command will issue an error message
    % if you try to execute this script on a computer without OpenGL.
    AssertOpenGL;

    % - Color Definitions: TODO: Use this code?
    %     black = BlackIndex(p.scrID);
    %     white = WhiteIndex(p.scrID);
    gray = 0.5;

    % pixels per degree
    p.fi = Screen('GetFlipInterval',w);   % fixation interval
    nf = round(p.stimDur/p.fi); % number of frames for gabor and noise
    p.stimDur = nf * p.fi; % update the real duration for record

    %Linearize Gamma:
    if p.linearize
        fid = fopen(p.gammaTableFile,'r');
        screen_clut = fread(fid,[256 3],'float64');
        fclose(fid);
        screen_clut = screen_clut - 1;
        screen_clut = screen_clut/255;
        Screen('LoadNormalizedGammaTable',p.scrID,screen_clut);
    else
        Screen('ReadNormalizedGammaTable',p.scrID);
    end

    HideCursor;
    ListenChar(2); % blocks input to matlab while function is running

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % generate fixation
    fix = ones(p.fixsz)*p.background/256;
    fix(p.fixsz/2 - p.fixwd/2 + 1: p.fixsz/2 + p.fixwd/2, :) = 0;
    fix(:, p.fixsz/2 - p.fixwd/2 + 1: p.fixsz/2 + p.fixwd/2) = 0;
    fix=Screen('MakeTexture',w,fix,[],[],2);

    %Instructions:
    Screen(w,'TextSize',28); Screen(w,'TextFont','Arial');
    DrawFormattedText(w, 'Indicate the direction of the gabor patches.\n\n\nPress the up arrow for up-tilted and the down arrow for down-tilted.\n\n\nPress any key to start.', 'center','center', 0);
    Screen(w,'Flip'); % show instruction

    %getResponse(space); % wait till any key to start
    KbPressWait(-3);
    Screen('Flip',w); % remove instruction

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Running the Experiment:
    p.start = datestr(now); % for record
    Priority(MaxPriority(w));

    %Create 'rec' Variable:
    % save trial number, location, and orientation in columns 1-3
    rec = nan(p.nTrials,8);
    rec(:,1) = (1:p.nTrials)';
    rec(:,7) = 99;

    for b = 1:p.nBlocks

        % Create a trial order with random order of trial types. Multiply matrix if
        % needed in units of trialTypeSize.
        trials = repmat(trialType,p.trialsPerBlock/trialTypeSize,1);
        trials = trials(randperm(size(trials,1)),:);
        rec(p.trialsPerBlock*(b-1)+1:p.trialsPerBlock*b,2:5) = trials;
        trialRange = p.trialsPerBlock * (b-1) + 1 : p.trialsPerBlock * b;
    
        for i = trialRange
            if rec(i,4) == -1 %stim loc, left or right
                stimrect = p.offsetRect(2, :); %stim on the right
                noiserect = p.offsetRect(1, :); %noise on the left
            elseif rec(i,4) == 1 %stim loc, left or right
                noiserect = p.offsetRect(2, :); %noise on the right
                stimrect = p.offsetRect(1, :); %stim on the left
            elseif rec(i,4) == 0 %noise noise trials
                stimrect = p.offsetRect(2, :); %stim on the right
                noiserect = p.offsetRect(1, :); %noise on the left
            end

            % generate gabor
            [x,y] = meshgrid(linspace(-p.radius, p.radius, p.halfPicSize * 2));
            gaborResult = exp(-((x.^2+y.^2)/2/p.sigma.^2)) .* ...
                sind(360* p.sf *(x*cosd(0) + y*sind(0))); % generates gabor matrix

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
            lowcutoff = 0.05; %TODO: add noise freq to parameters
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

            %noiseResult = Ntmp;
            noiseResult = envelope.*Ntmp;

            %Frames
            %slength = length(gaborResult);
            framexy = [-p.slength/2 p.slength/2 -p.slength/2 p.slength/2 ...
                -p.slength/2 -p.slength/2 p.slength/2 p.slength/2;...
                -p.slength/2 -p.slength/2 p.slength/2 p.slength/2 ...
                -p.slength/2 p.slength/2 -p.slength/2 p.slength/2];

            % Present inital fixation cross
            Screen('DrawTexture',w,fix,[],p.fixRect,0,1); % display fix
            %Screen('DrawLines',w,framexy,1,0,p.srectc(1,:),0); %left frame
            %Screen('DrawLines',w,framexy,1,0,p.srectc(2,:),0); %right frame
            Screen('Flip',w);
            WaitSecs(p.iti);

            % Present audio cue
            if rec(i,3) == 1
                sound(p.leftCue, 44100);
            elseif rec(i,3) == -1
                sound(p.rightCue, 44100);
            end

            WaitSecs(p.soa); %+ .001*randi(400)); % add random wait

            % create image matrixes
            gabor = gray + (gaborResult*p.gaborContrast + noiseResult*p.noiseContrast * gray);
            noiseOnly = gray + (noiseResult*p.noiseContrast + noiseResult*p.gaborContrast * gray);
            
            noiseTexture = Screen('MakeTexture',w,noiseOnly,[],[],2);
            stimTexture = Screen('MakeTexture',w,gabor,[],[],2);

            if rec(i, 4) == 0
                Screen('DrawTexture',w,noiseTexture,[],stimrect,p.dAngle*rec(i,5),1);
                Screen('DrawTexture',w,noiseTexture,[],noiserect,p.dAngle*rec(i,5),1);
            else
                Screen('DrawTexture',w,stimTexture,[],stimrect,p.dAngle*rec(i,5),1);
                Screen('DrawTexture',w,noiseTexture,[],noiserect,p.dAngle*rec(i,5),1);
            end
            Screen('DrawTexture',w,fix,[],p.fixRect,0,1); % display fix
            %Screen('DrawLines',w,framexy,1,0,p.srectc(1,:),0); %left frame
            %Screen('DrawLines',w,framexy,1,0,p.srectc(2,:),0); %right frame
            Screen('Flip',w);
            WaitSecs(p.stimDur);

            % start time and get response
            response = 0;
            tic

            Screen('DrawTexture',w,fix,[],p.fixRect,0,1); % display fix
            %Screen('DrawLines',w,framexy,1,0,p.srectc(1,:),0); %left frame
            %Screen('DrawLines',w,framexy,1,0,p.srectc(2,:),0); %right frame
            Screen('Flip',w);
           
            while toc < p.respDur;
                [down time kcode] = KbCheck(-3);

                % check for responses after stim presentation, record in rec
                if down
                    KeyPress = KbName(kcode);
                    if strcmp(KeyPress,'UpArrow')
                        response = 1;
                        rec(i,6) = response;
                        rec(i,7) = (response == rec(i,5));
                        rec(i,8) = toc;
                    elseif strcmp(KeyPress,'DownArrow')
                        response = -1;
                        rec(i,6) = response;
                        rec(i,7) = (response == rec(i,5));
                        rec(i,8) = toc;
                    elseif strcmp(KeyPress,'q')
                        error('esc:exit','ESC pressed.');
                    elseif strcmp(KeyPress, ~'DownArrow') || strcmp(KeyPress, ~'UpArrow') || strcmp(KeyPress, ~'q')
                        down = 0;
                    end
                end
%                 
                if rec(i,2) > 8 % if noise, response is false alarm
                    rec(i,7) = (down == false);
                    if abs(rec(i,6)) == 1
                        rec(i,7) = 0;
                    end
                end
            end
            
            % Feedback sounds
            if rec(i,7) == 1
                sound(p.beep);
            elseif rec(i,7) == 0 
                sound(p.beepError);
            else
                sound(p.beepError);
            end

            filename = sprintf('%s.mat', p.fname);
            cd Left_Right_Data
            save(filename,'p','rec');
            cd ..

        end % each trial, short analyzis
        
        if b < p.nBlocks
            t = 40*(1:p.nBlocks);
            block = rec(t-39:t,1:8);
            valid_mask = 5 <= block(:,2) & block(:,2) <= 8;
            valid_correct = length(block(valid_mask & block(:,7) == 1, 7));
            valid_incorrect = length(block(valid_mask & block(:,7) == 0, 7));
            valid_mean = 100 * valid_correct/(valid_correct+valid_incorrect);
            valid_rt = 1000 * mean(block(valid_mask & block(:,7) == 1, 8));
            
            pCorr = sprintf('%1.2g', valid_mean);
            rt = sprintf('%3.f', valid_rt);
            for t = 1:30
                DrawFormattedText(w, sprintf('Your score was %s%%.\n\nYour reaction time was %s ms. \n\nYou have completed %d out of %d blocks. Please take a break.\n\n\n%d', pCorr, rt, b, p.nBlocks, 31-t), 'center','center', 0);
                Screen('Flip', w);
                WaitSecs(1);
            end
            DrawFormattedText(w, sprintf('Your score was %s%%.\n\nYour reaction time was %s ms. \n\nYou have completed %d out of %d blocks. Please take a break.\n\n\nPress any key to continue.', pCorr, rt, b, p.nBlocks), 'center', 'center', 0);
            Screen('Flip', w);
            KbPressWait(-3);
        end
    end % each block
    
    cd Left_Right_Data
    ts = DATESTR(now,'yymmdd-HHMMSS');
    p.fnamefull = sprintf('%s_%s.mat', p.fname, ts);
    movefile(filename,p.fnamefull);
   % Run Stats 
   [s] = analyze_rec(rec,p.fnamefull);
   save(p.fnamefull,'p','rec','s');
    cd ..

    DrawFormattedText(w, 'Session Complete. Thank you!', 'center','center', 0);
    Screen('Flip', w);
    WaitSecs(2);

    Screen('CloseAll');
    Priority(0);

    ListenChar(0);

catch whatiserror % in case of error or user exit
    ShowCursor;
    Priority(0);
    ListenChar(0);
    Screen('Closeall');
    whatiserror.message
    whatiserror.stack
end

end % main function ends
