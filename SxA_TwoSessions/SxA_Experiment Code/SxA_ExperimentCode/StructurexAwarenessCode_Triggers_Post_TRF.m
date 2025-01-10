%% Structure x Awareness (Proj 1) Two Sessions
% Two session version of Structure and Awareness Project. First
% behavioural, second EEG. Includes Triggers. Needs
% instructions_session1.mat and instructions_session2.mat, as well as
% rhy.img and int.img in same folder.
%% Set-Up

% Clear
sca;
close all;
clear;
clear mex;

% Create Structs
screeninfo=[];
visStiminfo=[];
timeinfo=[];
textinfo=[];
mixinfo=[];


%% Main Parameters

experimentrun = str2num(cell2mat(inputdlg('Experiment Run? (1/0):'))); % adjusts all settings to run a whole experiment (including synctest, all practices etc.)
session_n = str2num(cell2mat(inputdlg('Session? (1/2):'))); 

conditions = [3 2 1]; % 1: rhythm, 2: interval, 3: irregular (later treated as blocks)
visStiminfo.maskintensity= 0.6; % default mask intensity, can be adjusted during threshold estimation
visStiminfo.gaborpercent=1-visStiminfo.maskintensity; % initial gabor intensity (will be adjusted during threshold)
visStiminfo.pracmaskintensity=0.4; % Very easy mask in the first ever practice round to make sure everyone knows what they are looking for
visStiminfo.pracgaborpercent=0.6;
mixinfo.startblock=1; % always start with block 1, unless you aborted midway, then this allows you to start from the block you want
timeinfo.maxresptime=3; % in seconds, after exceeding this response time, a warning will be shown

if ~experimentrun
    floatwin=0; % opens PTB in a floating window as opposed to taking over the entire screen
    screeninfo.scopetest=0; % disables also skip sync, speedrun, practice and threshold
    screeninfo.speedrun=0; % disables practice, threshold and quickpractice
    SkipSync=1;
    timeinfo.recalcframes=0; % recalculate frames each trial with current frame rate?
    PAScheck=1; %ask for reason why they chose that PAS response during practice

    % Rough Threshold
    visStiminfo.roughthreshold = 0; % run threshold estimation or run with standard mask contrast?
    thresholdsteps = 0.1;

    % Initial Practice
    initialpractice1=0; % objective practice
    initialpractice2=0; % adding subjective
    initialpractice3=0; % post-threshold practice with irregular pattern

    % Reminder Practice for 2nd Session
    reminderpractice=0;

    % 2AFC
    twoafc=0;

    % Catch trials for EEG?
    if session_n==1
    mixinfo.nCatchSample = 0; % catch trials per nContrasts (if 3 --> 36 catch trials for 12 samples)
    mixinfo.nCatchSampleirr = 0; % more samples (if 2 --> 40 catch trials for 12 samples)
    else
    mixinfo.nCatchSample = 3; % catch trials per nContrasts (if 3 --> 36 catch trials for 12 samples)
    mixinfo.nCatchSampleirr = 2; % more samples (if 2 --> 40 catch trials for 12 samples)
    end

elseif experimentrun && session_n==1
    screeninfo.scopetest=0; % disables also skip sync, speedrun, practice and threshold
    screeninfo.speedrun=0; % disables practice, threshold and quickpractice
    SkipSync=0;
    timeinfo.recalcframes=1; % recalculate frames each trial with current frame rate?
    PAScheck=1; % ask for reason why they chose that PAS response during practice
    floatwin=0; % opens PTB in a floating window as opposed to taking over the entire screen

    % Rough Threshold
    visStiminfo.roughthreshold = 1; % run threshold estimation or run with standard mask contrast?
    thresholdsteps = 0.1;

    % Initial Practice
    initialpractice1=1; % practice with visible targets and irregular pattern
    initialpractice2=1; % practice with threshold targets and irregular pattern
    initialpractice3=1;

    % 2AFC
    twoafc=1;

    % Catch trials for EEG?
    mixinfo.nCatchSample = 0; % catch trials per nContrasts (if 3 --> 36 catch trials for 12 samples)
    mixinfo.nCatchSampleirr = 0; % more samples (if 2 --> 40 catch trials for 12 samples)
% 
    HideCursor;
elseif experimentrun && session_n==2
    screeninfo.scopetest=0; % disables also skip sync, speedrun, practice and threshold
    screeninfo.speedrun=0; % disables practice, threshold and quickpractice
    SkipSync=0;
    timeinfo.recalcframes=1; % recalculate frames each trial with current frame rate?
    PAScheck=0; %ask for reason why they chose that PAS response during practice
    floatwin=0; % opens PTB in a floating window as opposed to taking over the entire screen

    % Rough Threshold
    visStiminfo.roughthreshold = 1; % run threshold estimation or run with standard mask contrast?
    thresholdsteps = 0.1;

    % Initial Practice
    initialpractice1=0; % practicee with visible targets and irregular pattern
    initialpractice2=0; % practice with threshold targets and irregular pattern
    initialpractice3=1; % Practice with circles and varying visibility

    % 2AFC
    twoafc=0;

    % Catch trials for EEG?
    mixinfo.nCatchSample = 3; % catch trials per nContrasts (if 3 --> 36 catch trials for 12 samples)
    mixinfo.nCatchSampleirr = 2; % more samples (if 2 --> 40 catch trials for 12 samples)

   HideCursor;
end

% Initial Practice
PracticeContrast=10;% fully visible
maxpractrials=10; % after this, the practice ends, but can be repeated
minvistrials=1; % after this, the experimenter can decide to skip the rest
minvistrials3=5; % for practice 3, after this, the experimenter can decide to skip the rest

% Practice before each block
nqp=2; % N quick-practice trials before each block, has to be divisible by 2

% Trial matrix function input
mixinfo.nContrasts = 10; % Number of contrasts excluding catch trial contrast (0), has to be 10, leave it like this

mixinfo.nSamples = 12; % how often each point is sampled (needs to be divisible by 2, currently doesn't affect trial matrix)
mixinfo.maxconscatch=1; % maximum consecutive catch trials
mixinfo.maxorientrep=10; % maximum trials of same orientation in a row (basically deactivated for now)
mixinfo.maxshortrep=2; % only irregular: maximum trials of short target time bin in a row
mixinfo.maxlongrep=2; % only irregular: maximum trials of long target time bin in a row
twoafcinfo.afcconseqcorr=4; %How many afc trials need to be correct consecutively
twoafcinfo.maxafcrep=15; % Maxmimum AFC trials
twoafcinfo.minafceval=4; % Minmum AFC trials, needs to be equal or more than afcconseqcorr


% Check trial split
if mod((mixinfo.nContrasts+mixinfo.nCatchSample)*mixinfo.nSamples,4)~=0
    error('Invalid combination of parameters.')
else
end

% Disable skip synch for scope test
if screeninfo.scopetest
    SkipSync=0;
    floatwin=0;
end
%% Secondary Parameters
% Timing in Seconds
% if ~screeninfo.speedrun
timeinfo.rhyISI = 0.8; % rhythm interval length
timeinfo.intISI = 0.8; % interval length
% else
%     rhyISI=0.1;
%     intISI=0.1;
% end

% Load externally stored instructions for current session (not in script to
% save space)
if session_n==1
    load('instructions_session1.mat')
else
    load('instructions_session2.mat')
end

irrmax=1.2;
irrmin=0.4;
timeinfo.irrwin = irrmin:0.05:irrmax; % irregular min and max
timeinfo.mindiff= 0.2; % minimum difference between irregular interval lengths (in seconds)
timeinfo.CueDur = 0.1; % Duration of cue
timeinfo.TarDur = 0.022; % Duration of target
timeinfo.TargetQuestInt = 1; % Interval between target and question
timeinfo.ITImask = [1:0.1:1.6];
timeinfo.prethr = [0.75 1 1.25]; % time before target during thresholding in seconds
timeinfo.preqjitter = [0.85 1 1.15]; % time after target before first question
%timeinfo.irrtartime = [timeinfo.rhyISI*0.8 timeinfo.rhyISI*0.9 timeinfo.rhyISI timeinfo.rhyISI*1.1 timeinfo.rhyISI*1.2]; % Irregular target times
timeinfo.irrtartime = [0.4 0.6 0.8 1 1.2]; % new version with increased jitter
timeinfo.ITI=1.4; % intertrial interval

% Calculations
mixinfo.nBlocks=3;
mixinfo.nBlocksirr=4;
mixinfo.nBlockstot=mixinfo.nBlocks*((length(conditions))-1)+mixinfo.nBlocksirr;
mixinfo.regularblength=((mixinfo.nContrasts+mixinfo.nCatchSample)*mixinfo.nSamples)/mixinfo.nBlocks;
mixinfo.irregularblength=((((mixinfo.nContrasts+mixinfo.nCatchSampleirr)*mixinfo.nSamples)/3)*5)/mixinfo.nBlocksirr; % 120 trials split across three time points, with two buffer time points on each side (=5)

% Gabor Information
sigma = 0.4;
numCycles = 21; % keep an odd number?
degrees=pi/180; % units
theta=45*degrees; % tilt

% Noise Parameters
visStiminfo.rectSize=300;
visStiminfo.noiseMean= 50;

% Response Prompts
textinfo.VisRespPrompt = 'Visibility?\n \n  (0) No experience\n (1) Brief glimpse\n (2) Almost clear experience\n (3) Clear experience';
textinfo.PAS1='0';
textinfo.PAS2='1';
textinfo.PAS3='2';
textinfo.PAS4='3';

% 2AFC Parameters
twoafcinfo.correct_timing=0.8;
twoafcinfo.offbeat_timing=[0.8*0.5 0.8*2]; % too short and too long
twoafcinfo.corrint=[1 2]; % Correct in 1st or second option

textinfo.warning=sprintf('Please remember to answer within  %.0f seconds after the question is shown.',timeinfo.maxresptime);

%% Subject Prompt and PTB Set Up

dirContent=dir;
subInFolder=[];
if ~screeninfo.speedrun
    subj = str2num(cell2mat(inputdlg('Enter Subject Number:')));
else
    subj=randi(100);
end

possiblesubj= 1:300;
if ~ismember(subj,possiblesubj)
    error('Invalid Subject Number')
elseif experimentrun
    fileNameAll=sprintf('SxA1.1_s%i-session%i',subj,session_n);
    for i=3:size(dirContent,1)
        currentFile = dirContent(i).name;
        subInFolder=[subInFolder, strcmp(currentFile  ,[fileNameAll '.mat'])];
    end

    if sum(subInFolder)>0
        error('subject num already exists')
    else
    end
elseif ~experimentrun
    fileNameAll=sprintf('test_%s', datestr(now,'mm-dd-yyyy_HH-MM'));
    for i=3:size(dirContent,1)
        currentFile = dirContent(i).name;
        subInFolder=[subInFolder, strcmp(currentFile  ,[fileNameAll '.mat'])];
    end
end

if experimentrun && session_n==1
    demogr.age = str2num(cell2mat(inputdlg('Age:'))); 
    demogr.hand = str2num(cell2mat(inputdlg('Handedness? (r/l):')));
end
demogr.subj=subj;
save([fileNameAll '.mat'], 'demogr')

% Trigger Set Up

timeinfo.triggerPortAddress=hex2dec('FFF8');

timeinfo.triggerPort=io64;
s=io64(timeinfo.triggerPort)
pause(1)
io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, 0);

% Setup PTB
PsychDefaultSetup(2);
Priority(1);
Screen('Preference', 'SkipSyncTests', SkipSync);
screeninfo.number = max(Screen('Screens'));

% Screen Values
white = WhiteIndex(screeninfo.number);
screeninfo.black = BlackIndex(screeninfo.number);
screeninfo.grey = white / 2;

if ~screeninfo.scopetest
    background=screeninfo.grey;
    screeninfo.fontcolour=screeninfo.black;
    if screeninfo.speedrun
        visStiminfo.roughthreshold = 0; % threshold is not coded to handle a speedrun
        initialpractice1=0; % practice with visible targets and irregular pattern
        initialpractice2=0; % practice with threshold targets and irregular pattern
        initialpractice3=0;
        twoafc=0;
        nqp=0;
    end
else
    background=screeninfo.black;
    screeninfo.fontcolour=white;
    visStiminfo.maskintensity= 0.01;
    visStiminfo.roughthreshold = 0; % run threshold estimation or run with standard mask contrast?
    initialpractice1=0; % practice with visible targets and irregular pattern
    initialpractice2=0; % practice with threshold targets and irregular pattern
    initialpractice3=0;
    twoafc=0;
    screeninfo.speedrun=0;
    nqp=0;
end     

if ~floatwin
[screeninfo.win, screeninfo.windowRect] = PsychImaging('OpenWindow', screeninfo.number, background, [], 32, 2,...
    [], [],  kPsychNeed32BPCFloat);
else
    [screeninfo.win, screeninfo.windowRect] = PsychImaging('OpenWindow', screeninfo.number, background,...
    [100 100 1600 1600], [], [], [], [], [], kPsychGUIWindow);
end

[screeninfo.xCenter, screeninfo.yCenter] = RectCenter(screeninfo.windowRect);
[screeninfo.axisx]=screeninfo.windowRect([1,3]);
[screeninfo.axisy]=screeninfo.windowRect([2,4]);
visStiminfo.baseRect = [0 0 120 120];
Screen('BlendFunction', screeninfo.win, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

Screen('TextFont', screeninfo.win, 'Ariel');
Screen('TextSize', screeninfo.win, 15);

screeninfo.ifi = Screen('GetFlipInterval', screeninfo.win);

% % Check frame rate
% curfrate=Screen('NominalFrameRate', screeninfo.win);
% if curfrate~=60
% %     error('Incorrect frame rate.')
%                                 sca
%                                 return
% end

%set(0, 'DefaultFigurePosition', [430 230 560 420])
%% Stimulus Information

% Fixation Cross
fixCrossDimPix = 10;
xCoords = [-fixCrossDimPix fixCrossDimPix 0 0];
yCoords = [0 0 -fixCrossDimPix fixCrossDimPix];
visStiminfo.fixCoords = [xCoords; yCoords];
visStiminfo.lineWidthPix = 2;

if screeninfo.scopetest
    visStiminfo.rectColorCue = 1;
end

% % Load Images
[rhyimg, rhyColorMap] = imread("rhyimg.png");
[intimg, intColorMap] = imread("intimg.png");

% Turn into textures
rhyimgtex = Screen('MakeTexture', screeninfo.win, rhyimg);
intimgtex = Screen('MakeTexture', screeninfo.win, intimg);
%% Experiment Prep
% Make condition matrix and counterbalance
conditionmatrix=[];
flipcondition=conditions;
flipcondition([2 3])=flipcondition([3 2]);

for c=1:mixinfo.nBlocks
    if mod(c,2)>0
        conditionmatrix=[conditionmatrix; conditions'];
    else
        conditionmatrix=[conditionmatrix; flipcondition'];
    end
end

conditionmatrix=[conditionmatrix; 3]; % add one more irregular block at the end (hard coded here unfortunately)

% Counterbalance
if mod(subj,2)>0
    conditionmatrix=changem(conditionmatrix,0,2);
    conditionmatrix=changem(conditionmatrix,2,1);
    conditionmatrix=changem(conditionmatrix,1,0);
else
end

% if mod(subj,4)==0 || mod(subj,4)==1
textinfo.KeyLeft = 'l';
textinfo.KeyRight = 'r';
textinfo.OrientRespPrompt ='Left or Right?';

% elseif mod(subj,4)==2 || mod(subj,4)==3
%     textinfo.KeyRight = 'a';
%     textinfo.KeyLeft = 'd';
%     textinfo.OrientRespPrompt ='Left or Right? \n \n  Right (a)        Left (d)';
%
% end


% Timing in Frames
timeinfo.rhyISIfr=round(timeinfo.rhyISI/screeninfo.ifi);
timeinfo.intISIfr=round(timeinfo.intISI/screeninfo.ifi);
timeinfo.CueDurFr=round(timeinfo.CueDur/screeninfo.ifi);
timeinfo.TarDurFr=round(timeinfo.TarDur/screeninfo.ifi);
timeinfo.TargQFr=round(timeinfo.TargetQuestInt/screeninfo.ifi);
timeinfo.ITIFr=round(timeinfo.ITImask/screeninfo.ifi);
timeinfo.thrFr=round(timeinfo.prethr/screeninfo.ifi);
timeinfo.irrtartimefr=round(timeinfo.irrtartime/screeninfo.ifi);

% 2afc:
twoafcinfo.correct_timing_fr=round(twoafcinfo.correct_timing/screeninfo.ifi);
twoafcinfo.offbeat_timing_fr=round(twoafcinfo.offbeat_timing/screeninfo.ifi); % too short and too long


% Other
visStiminfo.gcontrasts = [0 linspace(0.1,1,mixinfo.nContrasts)]; % log spaced contrasts for main experiment
visStiminfo.thresholdcontrasts = flip(0.1:thresholdsteps:1); % which contrasts are shown during staircase

% % Trial info
% targettrials=1:(length(visStiminfo.gcontrasts)-1); % minus one because we remove the catch trial contrast
% gaborconditions={[targettrials] [0 1]}; % Gabor contrasts(0 are catch trials, the other ones increasing contrasts) // 0-1: gabor orientation (0 right / 1 left)

% Make noise textures and cue/warning with default threshold
maxframes=3000; 
[x, y]=meshgrid(-1:2/(visStiminfo.rectSize-1):1,-1:2/(visStiminfo.rectSize-1):1);
wsLum=1; % luminance warning signal
csLum=0; % luminance cue signal
circleweak= 0.5; % how weak is the circle compared to mask (lower means less contrast)
gausscw = exp(-1.*(x.^2+y.^2)./sigma.^2); % gaussian envelope for cue/warning
gausscw=gausscw>0.3; % noisy but well defined circle (make bigger/smaller by changing this value)

DrawFormattedText(screeninfo.win, 'Loading. Please wait.', 'center', 'center', screeninfo.fontcolour);
Screen('Flip', screeninfo.win);

for fridx=1:maxframes
    if~screeninfo.scopetest
        visStiminfo.noiseimg=0.5+visStiminfo.maskintensity*(rand(visStiminfo.rectSize)-0.5);
        visStiminfo.pracnoiseimg=0.5+visStiminfo.pracmaskintensity*(rand(visStiminfo.rectSize)-0.5);
    else
        visStiminfo.noiseimg=0.01+visStiminfo.maskintensity*(rand(visStiminfo.rectSize)-0.5);
        visStiminfo.pracnoiseimg=0.01+visStiminfo.maskintensity*(rand(visStiminfo.rectSize)-0.5);
    end
    visStiminfo.noisetex(fridx)=Screen('MakeTexture', screeninfo.win, visStiminfo.noiseimg);
    visStiminfo.pracnoisetex(fridx)=Screen('MakeTexture', screeninfo.win, visStiminfo.pracnoiseimg);

    if ~screeninfo.scopetest
        % Cue and Warning (from NoisyCircleCode.m)
        randMaskW=(rand(visStiminfo.rectSize)-0.5)*2; %noisy mask
        randMaskC=(rand(visStiminfo.rectSize)-0.5)*2;

        rescaledWS=(2*wsLum-1)*gausscw;
        rescaledCS=(2*csLum-1)*gausscw;

        % blend outside
        warnStim=rescaledWS;
        cueStim=rescaledCS;

        warnStim(gausscw<0.3)=0.5+0.5*(visStiminfo.maskintensity*randMaskW(gausscw<0.3)+(1-visStiminfo.maskintensity)*rescaledWS(gausscw<0.3));
        cueStim(gausscw<0.3)=0.5+0.5*(visStiminfo.maskintensity*randMaskC(gausscw<0.3)+(1-visStiminfo.maskintensity)*rescaledCS(gausscw<0.3));

        % blend inside
        warnStim(gausscw>0.3)=0.5+0.5*(circleweak*randMaskW(gausscw>0.3)+(1-circleweak)*rescaledWS(gausscw>0.3));
        cueStim(gausscw>0.3)=0.5+0.5*(circleweak*randMaskC(gausscw>0.3)+(1-circleweak)*rescaledCS(gausscw>0.3));

    else

        % % Make Cues and Warning Signal Textures for scope test
        cueStim=0.01+visStiminfo.maskintensity*(rand(visStiminfo.rectSize)-0.5); %noisy mask
        xfrom=(length(cueStim)-visStiminfo.baseRect(3))/2;
        xto=xfrom+visStiminfo.baseRect(3);
        cueStim(xfrom:xto,xfrom:xto)=visStiminfo.rectColorCue;
        warnStim=cueStim; % all white for scope test

    end
    visStiminfo.warntex(fridx)=Screen('MakeTexture', screeninfo.win, warnStim);
    visStiminfo.cuetex(fridx)=Screen('MakeTexture', screeninfo.win, cueStim);
end

% Make Gabor Targets
[x, y]=meshgrid(-1:2/(visStiminfo.rectSize-1):1,-1:2/(visStiminfo.rectSize-1):1); %target

kxright=(pi*numCycles)*cos(theta);
kyright=(pi*numCycles)*sin(theta);
kxleft=(pi*numCycles)*cos(-theta);
kyleft=(pi*numCycles)*sin(-theta);

imRright=cos(kxright*x+kyright*y);
imRleft= cos(kxleft*x+kyleft*y);
imR = cos(2*numCycles*pi*(x)); %vertical

gaussEnvt = exp(-1.*(x.^2+y.^2)./sigma.^2); % create gaussian envelope (target)
gaussEnvt(gaussEnvt<0.01)=0; % change close-to-zero peripheral elements to zero
visStiminfo.gaborRight=imRright.*gaussEnvt;
visStiminfo.gaborLeft=imRleft.*gaussEnvt;
visStiminfo.gabor2afc=imR.*gaussEnvt; %vertical/neutral gabor for 2afc

% Create Gabor Example
exgableft=[screeninfo.axisx(2)*0.25, screeninfo.axisy(2)*0.75]; % loc of gabor
exgabright=[screeninfo.axisx(2)*0.75, screeninfo.axisy(2)*0.75]; % loc of gabor
gaborRect = [0 0 180 180];
gaborRectleft = CenterRectOnPointd(gaborRect, exgableft(1),exgableft(2));
gaborRectright = CenterRectOnPointd(gaborRect, exgabright(1),exgabright(2));

if ~screeninfo.scopetest
    exampletargetr = 0.5+0.5*(1*visStiminfo.gaborRight);
    exampletargetl = 0.5+0.5*(1*visStiminfo.gaborLeft);
else
    exampletargetr = 0.01+0.01*(1*visStiminfo.gaborRight);
    exampletargetl = 0.01+0.01*(1*visStiminfo.gaborLeft);
end

exampletargetr = Screen('MakeTexture', screeninfo.win, exampletargetr);
exampletargetl = Screen('MakeTexture', screeninfo.win, exampletargetl);

%% Text Locations

% Text Position of Instructions (not counterbalanced at the moment to avoid simon effect
% if mod(subj,4)==0 || mod(subj,4)==1
orientation1='Left';
orientation2='Right';
gaborright=gaborRectright;
gaborleft=gaborRectleft;

% Calculate Image and Text Positions
[vertsizerhy,~,~]=size(rhyimg);
[vertsizeint,~,~]=size(intimg);
rhytextpos=round((screeninfo.axisy(2)*0.6)+(vertsizerhy*0.5));
inttextpos=round((screeninfo.axisy(2)*0.6)+(vertsizeint*0.5));

[~,~,~,length1]=Screen('TextBounds', screeninfo.win, orientation1); % get bounds
[~,~,~,length2]=Screen('TextBounds', screeninfo.win, orientation2); % get bounds
textleft=[exgableft(1) round(exgableft(2)+gaborRect(4)*0.6)]; % shift on y axis in relation to gabor
textright=[exgabright(1) round(exgabright(2)+gaborRect(4)*0.6)]; % shift on y axis in relation to gabor
textleft(1)= textleft(1)-round(length1*0.5); % Centre text on location on x axis
textright(1)= textright(1)-round(length2*0.5); % Centre text on location on x axis

% Location for block instruction illustrations:
illustrationRect = [0 0 600 300];
illustrationRect = CenterRectOnPointd(illustrationRect, screeninfo.xCenter, screeninfo.yCenter+(0.25*screeninfo.axisy(2)));

%% Trial matrix
% Make one big trial matrix for the whole experiment

alltrials1=fulltrialmatrix(mixinfo,0); % rhythm trials
alltrials2=fulltrialmatrix(mixinfo,0); % interval trials
alltrials3=fulltrialmatrix(mixinfo,1); % irregular trials

% Save

subresults.trialmatrix1=alltrials1;
subresults.trialmatrix2=alltrials2;
subresults.trialmatrix3=alltrials3;
subresults.conditions=conditionmatrix;
save([fileNameAll '.mat'], 'subresults')
%% Run
try
    if ~screeninfo.scopetest
        %% Practice 1 (no circles, only orientation, fully visible, reduced mask contrast)
        while initialpractice1

            practicetrials1=[];

            % Practice Trial Matrix
            practicetrials1(1,1:(maxpractrials/2))=1;
            practicetrials1(1,(maxpractrials/2+1):maxpractrials)=0;
            [m,n] = size(practicetrials1);
            idx = randperm(n) ;
            b = practicetrials1 ;
            b(1,idx) = practicetrials1(1,:);
            practicetrials1=b';
            pcontrasts(1:maxpractrials)=PracticeContrast;
            allpracticetrials1=[pcontrasts' practicetrials1];

            % Present Instructions
            DrawFormattedText(screeninfo.win, textinfo.practice11, 'center', 'center',screeninfo.fontcolour);
            Screen('Flip', screeninfo.win);
            KbStrokeWait;
            DrawFormattedText(screeninfo.win, textinfo.practice12, 'center', 'center',screeninfo.fontcolour);
            Screen('DrawTexture', screeninfo.win, exampletargetr, [], gaborright);
            Screen('DrawTexture', screeninfo.win, exampletargetl, [], gaborleft);
            DrawFormattedText(screeninfo.win, orientation1, textleft(1),textleft(2),screeninfo.fontcolour);
            DrawFormattedText(screeninfo.win, orientation2, textright(1),textright(2),screeninfo.fontcolour);
            Screen('Flip', screeninfo.win);
            KbStrokeWait;
            %             textinfo.buttoninstr = sprintf('For left, press %s. \n\n For right. press %s. \n\n \n\n 3/4',textinfo.KeyLeft,textinfo.KeyRight);
            %             DrawFormattedText(screeninfo.win, textinfo.buttoninstr, 'center', 'center',screeninfo.fontcolour);
            %             Screen('Flip', screeninfo.win);
            %             KbStrokeWait;
            DrawFormattedText(screeninfo.win, textinfo.practice13, 'center', 'center',screeninfo.fontcolour);
            Screen('Flip', screeninfo.win);
            KbStrokeWait;

            praccorrect=0;
            contans=0;
            stoppr=0;

            % Run Very Visible Trials
            for nprac=1:maxpractrials
                % Run
                [~, ~, prOrRespEval, ~, ~, ~,~,~]=trialfunctionpro1(0,nprac,1,screeninfo, visStiminfo, timeinfo, textinfo, allpracticetrials1,0,0,1); % with reduced mask intensity
                % Give Feedback
                if prOrRespEval==1
                    DrawFormattedText(screeninfo.win,'Correct', 'center', 'center', screeninfo.fontcolour);
                    Screen('Flip', screeninfo.win);
                    praccorrect=praccorrect+1;
                    pause(1)
                elseif prOrRespEval==0
                    DrawFormattedText(screeninfo.win,'Incorrect', 'center', 'center', screeninfo.fontcolour);
                    Screen('Flip', screeninfo.win);
                    pause(1)
                end


                if nprac>=minvistrials
                    contans=0; % continue?-answer
                    while ~contans
                        DrawFormattedText(screeninfo.win, 'Keep practicing? (y-2/n-3/0-instr)', 'center', 'center',screeninfo.fontcolour);
                        Screen('Flip', screeninfo.win);
                        [~, keyNamethr, ~]=KbStrokeWait;
                        if strcmp(KbName(keyNamethr),'9')==1 || strcmp(KbName(keyNamethr),'9(')==1
                            DrawFormattedText(screeninfo.win,'Continue (8)? Exit (9)?', 'center', 'center', screeninfo.fontcolour);
                            Screen('Flip', screeninfo.win);
                            [~, keyNamethr, ~]=KbStrokeWait;
                            if strcmp(KbName(keyNamethr),'9') || strcmp(KbName(keyNamethr),'9(')==1
                                sca
                                return
                            elseif strcmp(KbName(keyNamethr),'8')|| strcmp(KbName(keyNamethr),'8*')==1
                                continue
                            end
                        elseif strcmp(KbName(keyNamethr), '3')==1 || strcmp(KbName(keyNamethr),'3#')==1
                            stoppr=1;
                            break
                        elseif strcmp(KbName(keyNamethr), '2')==1 || strcmp(KbName(keyNamethr),'2@')==1
                            contans=1;
                            %continue
                        elseif strcmp(KbName(keyNamethr),'0')==1
                            DrawFormattedText(screeninfo.win, textinfo.practice12, 'center', 'center',screeninfo.fontcolour);
                            Screen('DrawTexture', screeninfo.win, exampletargetr, [], gaborright);
                            Screen('DrawTexture', screeninfo.win, exampletargetl, [], gaborleft);
                            DrawFormattedText(screeninfo.win, orientation1, textleft(1),textleft(2),screeninfo.fontcolour);
                            DrawFormattedText(screeninfo.win, orientation2, textright(1),textright(2),screeninfo.fontcolour);
                            Screen('Flip', screeninfo.win);
                            KbStrokeWait;
                        end
                    end
                else
                end

                if stoppr
                    break
                end
            end

            praccorrect=0;
            contans=0;
            stoppr=0;

            DrawFormattedText(screeninfo.win, textinfo.practice14, 'center', 'center',screeninfo.fontcolour);
            Screen('Flip', screeninfo.win);
            KbStrokeWait;

            for nprac=1:maxpractrials
                % Run
                [~, ~, prOrRespEval, ~, ~, ~,~,~]=trialfunctionpro1(0,nprac,2,screeninfo, visStiminfo, timeinfo, textinfo, allpracticetrials1,0,0,1); % with reduced mask intensity
                % Give Feedback
                if prOrRespEval==1
                    DrawFormattedText(screeninfo.win,'Correct', 'center', 'center', screeninfo.fontcolour);
                    Screen('Flip', screeninfo.win);
                    praccorrect=praccorrect+1;
                    pause(1)
                elseif prOrRespEval==0
                    DrawFormattedText(screeninfo.win,'Incorrect', 'center', 'center', screeninfo.fontcolour);
                    Screen('Flip', screeninfo.win);
                    pause(1)
                end


                if nprac>=minvistrials
                    contans=0; % continue?-answer
                    while ~contans
                        DrawFormattedText(screeninfo.win, 'Keep practicing? (y-2/n-3/0-instr)', 'center', 'center',screeninfo.fontcolour);
                        Screen('Flip', screeninfo.win);
                        [~, keyNamethr, ~]=KbStrokeWait;
                        if strcmp(KbName(keyNamethr),'9')==1 || strcmp(KbName(keyNamethr),'9(')==1
                            DrawFormattedText(screeninfo.win,'Continue (8)? Exit (9)?', 'center', 'center', screeninfo.fontcolour);
                            Screen('Flip', screeninfo.win);
                            [~, keyNamethr, ~]=KbStrokeWait;
                            if strcmp(KbName(keyNamethr),'9') || strcmp(KbName(keyNamethr),'9(')==1
                                sca
                                return
                            elseif strcmp(KbName(keyNamethr),'8') || strcmp(KbName(keyNamethr),'8*')==1
                                continue
                            end
                        elseif strcmp(KbName(keyNamethr), '3')==1 || strcmp(KbName(keyNamethr),'3#')==1
                            stoppr=1;
                            break
                        elseif strcmp(KbName(keyNamethr), '2')==1 || strcmp(KbName(keyNamethr),'2@')==1
                            contans=1;
                            %continue
                        elseif strcmp(KbName(keyNamethr),'0')==1
                            DrawFormattedText(screeninfo.win, textinfo.practice12, 'center', 'center',screeninfo.fontcolour);
                            Screen('DrawTexture', screeninfo.win, exampletargetr, [], gaborright);
                            Screen('DrawTexture', screeninfo.win, exampletargetl, [], gaborleft);
                            DrawFormattedText(screeninfo.win, orientation1, textleft(1),textleft(2),screeninfo.fontcolour);
                            DrawFormattedText(screeninfo.win, orientation2, textright(1),textright(2),screeninfo.fontcolour);
                            Screen('Flip', screeninfo.win);
                            KbStrokeWait;
                        end
                    end
                else
                end

                if stoppr
                    break
                end
            end

            % Summary. Continue?
            while initialpractice1
                practicefeedback=sprintf('Practice 1 completed. \n Correct:%i/%i \n Repeat? y-2/n-3', praccorrect,nprac);
                DrawFormattedText(screeninfo.win, practicefeedback, 'center', 'center',screeninfo.fontcolour);
                Screen('Flip', screeninfo.win);
                [~, keyNamethr, ~]=KbStrokeWait;
                if strcmp(KbName(keyNamethr),'9')==1 || strcmp(KbName(keyNamethr),'9(')==1
                    DrawFormattedText(screeninfo.win,'Continue (8)? Exit (9)?', 'center', 'center', screeninfo.fontcolour);
                    Screen('Flip', screeninfo.win);
                    [~, keyNamethr, ~]=KbStrokeWait;
                    if strcmp(KbName(keyNamethr),'9') || strcmp(KbName(keyNamethr),'9(')==1
                        sca
                        return
                    elseif strcmp(KbName(keyNamethr),'8') || strcmp(KbName(keyNamethr),'8*')==1
                        continue
                    end
                elseif strcmp(KbName(keyNamethr), '3')==1 || strcmp(KbName(keyNamethr),'3#')==1
                    initialpractice1=0;
                    Screen('Flip', screeninfo.win);
                elseif strcmp(KbName(keyNamethr), '2')==1 || strcmp(KbName(keyNamethr),'2@')==1
                    break
                else
                end
            end
        end

        %% Practice 2 (no circles, both orientation and visibility, fully visible)
        while initialpractice2

            practicetrials2=[];

            % Practice Trial Matrix
            practicetrials2(1,1:(maxpractrials/2))=1;
            practicetrials2(1,(maxpractrials/2+1):maxpractrials)=0;
            [m,n] = size(practicetrials2);
            idx = randperm(n) ;
            b = practicetrials2 ;
            b(1,idx) = practicetrials2(1,:);
            practicetrials2=b';
            pcontrasts(1:maxpractrials)=PracticeContrast;
            allpracticetrials2=[pcontrasts' practicetrials2];

            % Present Instructions
            DrawFormattedText(screeninfo.win, textinfo.practice21, 'center', 'center',screeninfo.fontcolour);
            Screen('Flip', screeninfo.win);
            KbStrokeWait;
            DrawFormattedText(screeninfo.win, textinfo.practice22, 'center', 'center',screeninfo.fontcolour);
            Screen('Flip', screeninfo.win);
            KbStrokeWait;
            DrawFormattedText(screeninfo.win, textinfo.practice23, 'center', 'center',screeninfo.fontcolour);
            Screen('Flip', screeninfo.win);
            KbStrokeWait;

            praccorrect=0;
            contans=0;
            stoppr=0;

            % Run Visible Trials
            for nprac=1:maxpractrials
                % Run
                [~, ~, prOrRespEval, ~, ~, ~,~,~]=trialfunctionpro1(0,nprac,0,screeninfo, visStiminfo, timeinfo, textinfo, allpracticetrials2,0,0,1);
                % Give Feedback
                if prOrRespEval==1
                    DrawFormattedText(screeninfo.win,'Correct', 'center', 'center', screeninfo.fontcolour);
                    Screen('Flip', screeninfo.win);
                    praccorrect=praccorrect+1;
                    pause(1)
                elseif prOrRespEval==0
                    DrawFormattedText(screeninfo.win,'Incorrect', 'center', 'center', screeninfo.fontcolour);
                    Screen('Flip', screeninfo.win);
                    pause(1)
                end


                if nprac>=minvistrials
                    contans=0; % continue?-answer
                    while ~contans
                        DrawFormattedText(screeninfo.win, 'Keep practicing? (y-2/n-3)', 'center', 'center',screeninfo.fontcolour);
                        Screen('Flip', screeninfo.win);
                        [~, keyNamethr, ~]=KbStrokeWait;
                        if strcmp(KbName(keyNamethr),'9')==1 || strcmp(KbName(keyNamethr),'9(')==1
                            DrawFormattedText(screeninfo.win,'Continue (8)? Exit (9)?', 'center', 'center', screeninfo.fontcolour);
                            Screen('Flip', screeninfo.win);
                            [~, keyNamethr, ~]=KbStrokeWait;
                            if strcmp(KbName(keyNamethr),'9') || strcmp(KbName(keyNamethr),'9(')==1
                                sca
                                return
                            elseif strcmp(KbName(keyNamethr),'8') || strcmp(KbName(keyNamethr),'8*')==1
                                continue
                            end
                        elseif strcmp(KbName(keyNamethr), '3')==1 || strcmp(KbName(keyNamethr),'3#')==1
                            stoppr=1;
                            break
                        elseif strcmp(KbName(keyNamethr), '2')==1 || strcmp(KbName(keyNamethr),'2@')==1
                            contans=1;
                            %continue
                        else
                        end
                    end
                else
                end

                if stoppr
                    break
                end
            end

            % Summary. Continue?
            while initialpractice2
                practicefeedback=sprintf('Practice 1 completed. \n Correct:%i/%i \n Repeat? y-2/n-3', praccorrect,nprac);
                DrawFormattedText(screeninfo.win, practicefeedback, 'center', 'center',screeninfo.fontcolour);
                Screen('Flip', screeninfo.win);
                [~, keyNamethr, ~]=KbStrokeWait;
                if strcmp(KbName(keyNamethr),'9')==1 || strcmp(KbName(keyNamethr),'9(')==1
                    DrawFormattedText(screeninfo.win,'Continue (8)? Exit (9)?', 'center', 'center', screeninfo.fontcolour);
                    Screen('Flip', screeninfo.win);
                    [~, keyNamethr, ~]=KbStrokeWait;
                    if strcmp(KbName(keyNamethr),'9') || strcmp(KbName(keyNamethr),'9(')==1
                        sca
                        return
                    elseif strcmp(KbName(keyNamethr),'8') || strcmp(KbName(keyNamethr),'8*')==1
                        continue
                    end
                elseif strcmp(KbName(keyNamethr), '3')==1 || strcmp(KbName(keyNamethr),'3#')==1
                    initialpractice2=0;
                    Screen('Flip', screeninfo.win);
                elseif strcmp(KbName(keyNamethr), '2')==1 || strcmp(KbName(keyNamethr),'2@')==1
                    break
                else
                end
            end
        end
        
        %% Introduction to Session 2
        if session_n==2
                    DrawFormattedText(screeninfo.win,textinfo.introduction11, 'center', 'center', screeninfo.fontcolour);
                    Screen('Flip', screeninfo.win);
                    KbStrokeWait;

                    DrawFormattedText(screeninfo.win,textinfo.introduction12, 'center', 'center', screeninfo.fontcolour);
                    Screen('DrawTexture', screeninfo.win, exampletargetr, [], gaborright);
                    Screen('DrawTexture', screeninfo.win, exampletargetl, [], gaborleft);
                    DrawFormattedText(screeninfo.win, orientation1, textleft(1),textleft(2),screeninfo.fontcolour);
                    DrawFormattedText(screeninfo.win, orientation2, textright(1),textright(2),screeninfo.fontcolour);
                    Screen('Flip', screeninfo.win);
                    KbStrokeWait;

                    DrawFormattedText(screeninfo.win,textinfo.introduction13, 'center', 'center', screeninfo.fontcolour);
                    Screen('Flip', screeninfo.win);
                    KbStrokeWait;

                    DrawFormattedText(screeninfo.win,textinfo.introduction14, 'center', 'center', screeninfo.fontcolour);
                    Screen('Flip', screeninfo.win);
                    KbStrokeWait;
        end
        %% Threshold
        if visStiminfo.roughthreshold==1
            thresholdrevers=20;  % will crash if you put 1, so put at least 2, for experiment use 20
            threshrun=1;
            thrtrialoverall=1;
            roughthreshold=1;
            while roughthreshold==1
                thrtrial=1;
                thresholdfound=0;
                adjusted=0;
                thc=1; % counter
                prevstep=0; % up or down?
                currreverse=0; % how many reversals up until now?
                doublelow=0; % not seen twice at lowest contrast
                firstreversal=0; % at which trial did the first reversal happen?
                firstreveralrec=0; 
                if threshrun==1 %only show instructions for first threshold round
                    DrawFormattedText(screeninfo.win, textinfo.threshold1, 'center', 'center',screeninfo.fontcolour);
                    Screen('DrawTexture', screeninfo.win, exampletargetr, [], gaborright);
                    Screen('DrawTexture', screeninfo.win, exampletargetl, [], gaborleft);
                    DrawFormattedText(screeninfo.win, orientation1, textleft(1),textleft(2),screeninfo.fontcolour);
                    DrawFormattedText(screeninfo.win, orientation2, textright(1),textright(2),screeninfo.fontcolour);
                    Screen('Flip', screeninfo.win);
                    KbStrokeWait;
                else
                    DrawFormattedText(screeninfo.win, 'Press any button to start.', 'center', 'center',screeninfo.fontcolour);
                    Screen('Flip', screeninfo.win);
                    KbStrokeWait;
                end

                while ~thresholdfound
                    finishadjusting=0;
                    if thc<=length(visStiminfo.thresholdcontrasts) % as long as max trials has not been reached, loop over the contrasts
                        thrcurrcontrast= visStiminfo.thresholdcontrasts(thc);
                        [prOrRespEval, prVisResp,repeat]=thresholdfunctionpro1(thrcurrcontrast,screeninfo, visStiminfo, timeinfo, textinfo);

                        %                     if prOrRespEval==1
                        %                         DrawFormattedText(screeninfo.win,'Correct', 'center', 'center', screeninfo.fontcolour);
                        %                         Screen('Flip', screeninfo.win);
                        %                         %             thc=thc+1;
                        %                         pause(1)
                        %                     elseif prOrRespEval==0
                        %                         DrawFormattedText(screeninfo.win,'Incorrect', 'center', 'center', screeninfo.fontcolour);
                        %                         Screen('Flip', screeninfo.win);
                        %                         %             thc=thc+1;
                        %                         pause(1)
                        if repeat==1 % if threshold trial was interrupted, do not continue but repeat the trial instead
                            DrawFormattedText(screeninfo.win,'Invalid. Repeat.', 'center', 'center', screeninfo.fontcolour);
                            Screen('Flip', screeninfo.win);
                            pause(1)
                        end

                        % Record values
                        if currreverse==1 && ~firstreveralrec
                            firstreversal=thrtrial-1;
                            firstreveralrec=1;
                        end

                        subresults.thresholdresults(:,thrtrialoverall)=[threshrun; visStiminfo.gaborpercent; thrcurrcontrast;  prVisResp];
                        thrtrialoverall=thrtrialoverall+1;
                        thrtrial=thrtrial+1;

                        % Staircase (down in steps of 2 initially, steps of 1 after first reversal)
                        if currreverse==0 && ~prVisResp==0 && prOrRespEval && ~repeat % if first down, and response is seen and correct and trial was not interrupted
                            thc=thc+2; % go two contrasts lower
                            if thc>10 && doublelow==0
                                thc=10;
                                doublelow=1;
                            end
                            if strcmp(prevstep,'up') % if direction is reversed, add to counter
                                currreverse=currreverse+1;
                            end
                            prevstep='down';
                        elseif ~currreverse==0 && ~prVisResp==0 && prOrRespEval && ~repeat % if not first down, and response is seen and correct and not interrupted
                            thc=thc+1; % go one down
                            if thc>10 && doublelow==0
                                thc=10;
                                doublelow=1;
                            end
                            if strcmp(prevstep,'up') % if direction is reversed, add to counter
                                currreverse=currreverse+1;
                            end
                            prevstep='down';
                        elseif prVisResp==0 && ~repeat % if not seen and not interrupted
                            thc=thc-1; % go one up
                            if strcmp(prevstep,'down') % if direction is reversed, add to counter
                                currreverse=currreverse+1;
                            end
                            prevstep='up';
                        end

                        % Threshold reached, adjust gabor percentage
                        if prVisResp==0 && thc<1 % if participant does not see target at highest contrast
                            while ~finishadjusting==1
                                ImpossibleValueText=sprintf('Not seen at highest contrast. Adjust intensity? (y-2/n-3)');
                                DrawFormattedText(screeninfo.win, ImpossibleValueText, 'center', 'center',screeninfo.fontcolour);
                                Screen('Flip', screeninfo.win);
                                [~, keyNamethr, ~]=KbStrokeWait;
                                if strcmp(KbName(keyNamethr),'9')==1 || strcmp(KbName(keyNamethr),'9(')==1
                                    DrawFormattedText(screeninfo.win,'Continue (8)? Exit (9)?', 'center', 'center', screeninfo.fontcolour);
                                    Screen('Flip', screeninfo.win);
                                    [~, keyNamethr, ~]=KbStrokeWait;
                                    if strcmp(KbName(keyNamethr),'9') || strcmp(KbName(keyNamethr),'9(')==1
                                        sca
                                        return
                                    elseif strcmp(KbName(keyNamethr),'8') || strcmp(KbName(keyNamethr),'8*')==1
                                        continue
                                    end
                                elseif strcmp(KbName(keyNamethr), '2')==1 || strcmp(KbName(keyNamethr),'2@')==1
                                    Screen('Flip', screeninfo.win);
                                    recintensity=(visStiminfo.gaborpercent/10)*(1*10)*2;
                                    while ~finishadjusting
                                        contrasttext=sprintf('Current intensity: %.2f \n\n Recommended Intensity: %.2f \n\n Press any button to adjust.',visStiminfo.gaborpercent,recintensity);
                                        DrawFormattedText(screeninfo.win,contrasttext, 'center', 'center', screeninfo.fontcolour);
                                        Screen('Flip', screeninfo.win);
                                        KbStrokeWait;
                                        newcontrast=adjustintensity(screeninfo, visStiminfo.gaborpercent); % confirm new threshold
                                        Screen('Flip', screeninfo.win);
                                        Confirmationtext=sprintf('New intensity: %.2f \n Correct? (y-2/n-3)',newcontrast);
                                        DrawFormattedText(screeninfo.win,Confirmationtext, 'center', 'center', screeninfo.fontcolour);
                                        Screen('Flip', screeninfo.win);
                                        [~, keyNamethr, ~]=KbStrokeWait;
                                        if strcmp(KbName(keyNamethr),'9')==1 || strcmp(KbName(keyNamethr),'9(')==1
                                            DrawFormattedText(screeninfo.win,'Continue (8)? Exit (9)?', 'center', 'center', screeninfo.fontcolour);
                                            Screen('Flip', screeninfo.win);
                                            [~, keyNamethr, ~]=KbStrokeWait;
                                            if strcmp(KbName(keyNamethr),'9') || strcmp(KbName(keyNamethr),'9(')==1
                                                sca
                                                return
                                            elseif strcmp(KbName(keyNamethr),'8') || strcmp(KbName(keyNamethr),'8*')==1
                                                continue
                                            end
                                        elseif strcmp(KbName(keyNamethr), '2')==1 || strcmp(KbName(keyNamethr),'2@')==1 %continue
                                            finishadjusting=1;
                                            thresholdfound=1;
                                            visStiminfo.gaborpercent=newcontrast;
                                            adjusted=1;
                                        elseif strcmp(KbName(keyNamethr),'3')==1 || strcmp(KbName(keyNamethr),'3#')==1 % type threshold again
                                        end
                                    end
                                elseif strcmp(KbName(keyNamethr),'3')==1 || strcmp(KbName(keyNamethr),'3#')==1
                                    finishadjusting=1;
                                    thresholdfound=1;
                                else
                                end
                            end
                            break
                        elseif currreverse==thresholdrevers % if staircase has turned around N times
                            while ~finishadjusting==1
                                % Calculate threshold (average from first reversal to end)
                                % subresults.calculatedthreshold(threshrun)=mean(subresults.thresholdresults(3,end-6:end));
                                subresults.firstreversal(threshrun)=firstreversal;
                                resultscurrentrun=subresults.thresholdresults(3,subresults.thresholdresults(1,:)==threshrun);
                                subresults.calculatedthreshold(threshrun)=mean(resultscurrentrun(firstreversal:end));
                                Confirmationtext=sprintf('Threshold reached. \n\n Threshold: %.2f \n\n Adjust intensity? (y-2/n-3)',subresults.calculatedthreshold(threshrun));
                                DrawFormattedText(screeninfo.win,Confirmationtext, 'center', 'center',screeninfo.fontcolour);
                                Screen('Flip', screeninfo.win);
                                [~, keyNamethr, ~]=KbStrokeWait;
                                if strcmp(KbName(keyNamethr),'9')==1 || strcmp(KbName(keyNamethr),'9(')==1
                                    DrawFormattedText(screeninfo.win,'Continue (8)? Exit (9)?', 'center', 'center', screeninfo.fontcolour);
                                    Screen('Flip', screeninfo.win);
                                    [~, keyNamethr, ~]=KbStrokeWait;
                                    if strcmp(KbName(keyNamethr),'9') || strcmp(KbName(keyNamethr),'9(')==1
                                        sca
                                        return
                                    elseif strcmp(KbName(keyNamethr),'8') || strcmp(KbName(keyNamethr),'8*')==1
                                        continue
                                    end
                                elseif strcmp(KbName(keyNamethr), '2')==1 || strcmp(KbName(keyNamethr),'2@')==1
                                    Screen('Flip', screeninfo.win);
                                    recintensity=(visStiminfo.gaborpercent/10)*(subresults.calculatedthreshold(threshrun)*10)*2;
                                    while ~finishadjusting
                                        newcontrast=adjustintensity(screeninfo, visStiminfo.gaborpercent);
                                        Screen('Flip', screeninfo.win);
                                        Confirmationtext=sprintf('New intensity: %.2f \n\n Correct? (y-2/n-3)',newcontrast);
                                        DrawFormattedText(screeninfo.win,Confirmationtext, 'center', 'center', screeninfo.fontcolour);
                                        Screen('Flip', screeninfo.win);
                                        [~, keyNamethr, ~]=KbStrokeWait;
                                        if strcmp(KbName(keyNamethr),'9')==1 || strcmp(KbName(keyNamethr),'9(')==1
                                            DrawFormattedText(screeninfo.win,'Continue (8)? Exit (9)?', 'center', 'center', screeninfo.fontcolour);
                                            Screen('Flip', screeninfo.win);
                                            [~, keyNamethr, ~]=KbStrokeWait;
                                            if strcmp(KbName(keyNamethr),'9') || strcmp(KbName(keyNamethr),'9(')==1
                                                sca
                                                return
                                            elseif strcmp(KbName(keyNamethr),'8') || strcmp(KbName(keyNamethr),'8*')==1
                                                continue
                                            end
                                        elseif strcmp(KbName(keyNamethr), '2')==1 || strcmp(KbName(keyNamethr),'2@')==1 %continue
                                            visStiminfo.gaborpercent=newcontrast;
                                            finishadjusting=1;
                                            thresholdfound=1;
                                            adjusted=1;
                                        elseif strcmp(KbName(keyNamethr),'3')==1 || strcmp(KbName(keyNamethr),'3#')==1 % type threshold again
                                        end
                                    end
                                elseif strcmp(KbName(keyNamethr),'3')==1 || strcmp(KbName(keyNamethr),'3#')==1
                                    finishadjusting=1;
                                    thresholdfound=1;
                                else
                                end
                            end
                            break
                        elseif thc>length(visStiminfo.thresholdcontrasts) && doublelow % if lowest contrast has been detected two times
                            while ~finishadjusting==1
                                DrawFormattedText(screeninfo.win,'Detection at all contrasts. Adjust intensity? (y-2/n-3)', 'center', 'center', screeninfo.fontcolour);
                                Screen('Flip', screeninfo.win);
                                [~, keyNamethr, ~]=KbStrokeWait;
                                if strcmp(KbName(keyNamethr),'9')==1 || strcmp(KbName(keyNamethr),'9(')==1
                                    DrawFormattedText(screeninfo.win,'Continue (8)? Exit (9)?', 'center', 'center', screeninfo.fontcolour);
                                    Screen('Flip', screeninfo.win);
                                    [~, keyNamethr, ~]=KbStrokeWait;
                                    if strcmp(KbName(keyNamethr),'9') || strcmp(KbName(keyNamethr),'9(')==1
                                        sca
                                        return
                                    elseif strcmp(KbName(keyNamethr),'8') || strcmp(KbName(keyNamethr),'8*')==1
                                        continue
                                    end
                                elseif strcmp(KbName(keyNamethr), '2')==1 || strcmp(KbName(keyNamethr),'2@')==1 % enter new threshold
                                    recintensity=(visStiminfo.gaborpercent/10)*(1*10)*2;
                                    contrasttext=sprintf('Current intensity: %.2f \n\n Recomended Intensity: %.2f \n\n Press any button to adjust.',visStiminfo.gaborpercent,recintensity);
                                    DrawFormattedText(screeninfo.win,contrasttext, 'center', 'center', screeninfo.fontcolour);
                                    Screen('Flip', screeninfo.win);
                                    KbStrokeWait;
                                    newcontrast=adjustintensity(screeninfo, visStiminfo.gaborpercent);
                                    Screen('Flip', screeninfo.win);
                                    Confirmationtext=sprintf('New intensity: %.2f \n\n Correct? (y/n)',newcontrast);
                                    DrawFormattedText(screeninfo.win,Confirmationtext, 'center', 'center', screeninfo.fontcolour);
                                    Screen('Flip', screeninfo.win);
                                    [~, keyNamethr, ~]=KbStrokeWait;
                                    if strcmp(KbName(keyNamethr),'9')==1 || strcmp(KbName(keyNamethr),'9(')==1
                                        DrawFormattedText(screeninfo.win,'Continue (8)? Exit (9)?', 'center', 'center', screeninfo.fontcolour);
                                        Screen('Flip', screeninfo.win);
                                        [~, keyNamethr, ~]=KbStrokeWait;
                                        if strcmp(KbName(keyNamethr),'9') || strcmp(KbName(keyNamethr),'9(')==1
                                            sca 
                                            return
                                        elseif strcmp(KbName(keyNamethr),'3')
                                            continue
                                        end
                                    elseif strcmp(KbName(keyNamethr), '2')==1 || strcmp(KbName(keyNamethr),'2@')==1 %continue
                                        visStiminfo.gaborpercent=newcontrast;
                                        finishadjusting=1;
                                        thresholdfound=1;
                                        adjusted=1;
                                    elseif strcmp(KbName(keyNamethr),'3')==1 || strcmp(KbName(keyNamethr),'3#')==1 % type threshold again
                                    end
                                elseif strcmp(KbName(keyNamethr),'3')==1 || strcmp(KbName(keyNamethr),'3#')==1
                                    finishadjusting=1;
                                    thresholdfound=1;
                                    break
                                else
                                end
                            end
                        else
                        end
                    end
                end

                % One round of thresholding finished. Repeat
                answered=0;
                while ~answered
                    DrawFormattedText(screeninfo.win,'Repeat threshold? (y-2/n-3)', 'center', 'center',screeninfo.fontcolour);
                    Screen('Flip', screeninfo.win);
                    [~, keyNamethr, ~]=KbStrokeWait;
                    if strcmp(KbName(keyNamethr),'9')==1 || strcmp(KbName(keyNamethr),'9(')==1
                        DrawFormattedText(screeninfo.win,'Continue (8)? Exit (9)?', 'center', 'center', screeninfo.fontcolour);
                        Screen('Flip', screeninfo.win);
                        [~, keyNamethr, ~]=KbStrokeWait;
                        if strcmp(KbName(keyNamethr),'9') || strcmp(KbName(keyNamethr),'9(')==1
                            sca
                            return
                        elseif strcmp(KbName(keyNamethr),'8') || strcmp(KbName(keyNamethr),'8*')==1
                            continue
                        end
                    elseif strcmp(KbName(keyNamethr),'3')==1 || strcmp(KbName(keyNamethr),'3#')==1
                        roughthreshold=0;
                        answered=1;
                    elseif strcmp(KbName(keyNamethr), '2')==1 || strcmp(KbName(keyNamethr),'2@')==1
                        threshrun=threshrun+1;
                        answered=1;
                    end
                end
            end

            %         % Plot threshold values
            %         figure;
            %         for i=1:threshrun
            %             subplot(2,2,i)
            %             runidx=subresults.thresholdresults(1,:)==i;
            %             plot(subresults.thresholdresults(3,runidx),'-o')
            %             plottitle=sprintf('Gabor Percent: %.2f /// Threshold: %.2f',mean(subresults.thresholdresults(2,runidx)),subresults.calculatedthreshold(i));
            %             title(plottitle)
            %             ylim([0 1])
            %             yticks(flip(visStiminfo.thresholdcontrasts))
            %             hold on
            %             yline(subresults.calculatedthreshold(i),'r','LineWidth',2) % plot the threshold (average of the last three turns)
            %             hold off
            %         end
        end
        %% Practice 3 (including PAS check)

        while initialpractice3

            % Practice Trial Matrix
            practicetrials3=[];

            %create left and Right orientations
            practicetrials3(1,1:(maxpractrials/2))=1;
            practicetrials3(1,(maxpractrials/2+1):maxpractrials)=0;

            % Shuffle
            [m,n] = size(practicetrials3);
            idx = randperm(n);
            b = practicetrials3 ;
            b(1,idx) = practicetrials3(1,:);
            practicetrials3=b';

            % Add contrasts
            pcontrasts=1:mixinfo.nContrasts;
            allpracticetrials3=[pcontrasts' practicetrials3];

            % Shuffle again
            [m,n] = size(allpracticetrials3);
            idx = randperm(m)';
            b = allpracticetrials3 ;
            allpracticetrials3(idx,1) = allpracticetrials3(:,1);

            % Make sure that the first one is always full visibility (there
            % is probably a more elegant way to do this, but this works)
            idx_fullcont=find(allpracticetrials3(:,1)==10);
            allpracticetrials3 = [allpracticetrials3(idx_fullcont,:);allpracticetrials3];
            allpracticetrials3(idx_fullcont+1,:)=[];

            % Present Instructions
            if session_n==1
            DrawFormattedText(screeninfo.win, textinfo.practice31, 'center', 'center',screeninfo.fontcolour);
            elseif session_n==2
            DrawFormattedText(screeninfo.win, 'Let`s do a few practice trials. \n\n\n\n Press any button to start', 'center', 'center',screeninfo.fontcolour);
            end
            Screen('Flip', screeninfo.win);
            KbStrokeWait;

            praccorrect=0;

            % Run Threshold Practice Trials

            % PAS counters
            PAScount0=0;
            PAScount1=0;
            PAScount2=0;
            PAScount3=0;

            for nprac=1:maxpractrials
                % Run
                [~, ~, prOrRespEval, PASresp, ~, ~, ~,~]=trialfunctionpro1(0,nprac,0,screeninfo, visStiminfo, timeinfo, textinfo, allpracticetrials3,0,1,0);
                % Give Feedback
                if prOrRespEval==1
                    DrawFormattedText(screeninfo.win,'Correct', 'center', 'center', screeninfo.fontcolour);
                    Screen('Flip', screeninfo.win);
                    praccorrect=praccorrect+1;
                    pause(1)
                elseif prOrRespEval==0
                    DrawFormattedText(screeninfo.win,'Incorrect', 'center', 'center', screeninfo.fontcolour);
                    Screen('Flip', screeninfo.win);
                    pause(1)
                end

                if PASresp==0
                    currPAScount=PAScount0;
                    PAScount0=PAScount0+1;
                elseif PASresp==1
                    currPAScount=PAScount1;
                    PAScount1=PAScount1+1;
                elseif PASresp==2
                    currPAScount=PAScount2;
                    PAScount2=PAScount2+1;
                elseif PASresp==3
                    currPAScount=PAScount3;
                    PAScount3=PAScount3+1;
                end

                if PAScheck && currPAScount==1
                    if PASresp==0
                        DrawFormattedText(screeninfo.win,'For visibility you responded "No experience", please explain why you chose this option.', 'center', 'center', screeninfo.fontcolour);
                    elseif PASresp==1
                        DrawFormattedText(screeninfo.win,'For visibility you responded "Brief glimpse", please explain why you chose this option.', 'center', 'center', screeninfo.fontcolour);
                    elseif PASresp==2
                        DrawFormattedText(screeninfo.win,'For visibility you responded "Almost clear experience", please explain why you chose this option.', 'center', 'center', screeninfo.fontcolour);
                    elseif PASresp==3
                        DrawFormattedText(screeninfo.win,'For visibility you responded "Clear experience", please explain why you chose this option.', 'center', 'center', screeninfo.fontcolour);
                    end
                    Screen('Flip', screeninfo.win);
                    KbStrokeWait;
                    KbStrokeWait;
                    DrawFormattedText(screeninfo.win,'Press any key to start the next practice trial.', 'center', 'center', screeninfo.fontcolour);
                    Screen('Flip', screeninfo.win);
                    KbStrokeWait;
                end

 
%                 For a participant to be exposed to all contrasts, it has to
%                 run all the trials. if the experimenter should be able to
%                 skip part of the practice, activate the code below.
                            stoppr=0;
                            if nprac>=minvistrials3
                                contans=0;
                                while ~contans
                                DrawFormattedText(screeninfo.win, 'Keep practicing? (y2/n3)', 'center', 'center',screeninfo.fontcolour);
                                Screen('Flip', screeninfo.win);
                                [~, keyNamethr, ~]=KbStrokeWait;
                                if strcmp(KbName(keyNamethr),'e')==1
                                    DrawFormattedText(screeninfo.win,'Continue (c)? Exit (e)?', 'center', 'center', screeninfo.fontcolour);
                                    Screen('Flip', screeninfo.win);
                                    [~, keyNamethr, ~]=KbStrokeWait;
                                    if strcmp(KbName(keyNamethr),'e')
                                        sca
                                        return
                                    elseif strcmp(KbName(keyNamethr),'c')
                                        continue
                                    end
                                elseif strcmp(KbName(keyNamethr), '3')==1 || strcmp(KbName(keyNamethr),'3#')==1
                                    stoppr=1;
                                    break
                                elseif strcmp(KbName(keyNamethr), '2')==1 || strcmp(KbName(keyNamethr),'2@')==1
                                    contans=1;
                                    %continue
                                else
                                end
                                end
                            else
                            end

                            if stoppr
                                break
                            end
            end

            % Summary. Continue?
            while initialpractice3
                practicefeedback=sprintf('Practice completed. \n Correct:%i/%i \n Repeat? y-2/n-3', praccorrect,nprac);
                DrawFormattedText(screeninfo.win, practicefeedback, 'center', 'center',screeninfo.fontcolour);
                Screen('Flip', screeninfo.win);
                [~, keyNamethr, ~]=KbStrokeWait;
                if strcmp(KbName(keyNamethr),'9')==1 || strcmp(KbName(keyNamethr),'9(')==1
                    DrawFormattedText(screeninfo.win,'Continue (8)? Exit (9)?', 'center', 'center', screeninfo.fontcolour);
                    Screen('Flip', screeninfo.win);
                    [~, keyNamethr, ~]=KbStrokeWait;
                    if strcmp(KbName(keyNamethr),'9') || strcmp(KbName(keyNamethr),'9(')==1
                        sca
                        return
                    elseif strcmp(KbName(keyNamethr),'8') || strcmp(KbName(keyNamethr),'8*')==1
                        continue
                    end
                elseif strcmp(KbName(keyNamethr), '3')==1 || strcmp(KbName(keyNamethr),'3#')==1
                    initialpractice3=0;
                    Screen('Flip', screeninfo.win);
                elseif strcmp(KbName(keyNamethr), '2')==1 || strcmp(KbName(keyNamethr),'2@')==1
                    break
                else
                end
            end
        end
        %% Run Experiment

        %Instruction Screen
        DrawFormattedText(screeninfo.win, textinfo.instructions, 'center', 'center',screeninfo.fontcolour);
        Screen('Flip', screeninfo.win);
        KbStrokeWait;

        if session_n==2
        DrawFormattedText(screeninfo.win, textinfo.instructions2, 'center', 'center',screeninfo.fontcolour);
        Screen('Flip', screeninfo.win);
        KbStrokeWait;

        DrawFormattedText(screeninfo.win, textinfo.instructions3, 'center', 'center',screeninfo.fontcolour);
        Screen('Flip', screeninfo.win);
        KbStrokeWait;
        end

        subresults.visStimifo=visStiminfo;
        subresults.mixinfo=timeinfo;
        subresults.timeinfo=mixinfo;
 
%         startexp=0;

%         while ~startexp
%             DrawFormattedText(screeninfo.win, 'Start experiment? (y-2)\nShow instructions again? (i-5)', 'center', 'center',screeninfo.fontcolour);
%             Screen('Flip', screeninfo.win);
%             answered=0;
%             while ~answered
%                 [~, keyNamePrac, ~]=KbStrokeWait;
%                 if strcmp(KbName(keyNamePrac),'9')==1 || strcmp(KbName(keyNamethr),'9(')==1
%                     DrawFormattedText(screeninfo.win,'Continue (8)? Exit (9)?', 'center', 'center', screeninfo.fontcolour);
%                     Screen('Flip', screeninfo.win);
%                     [~, keyNamethr, ~]=KbStrokeWait;
%                     if strcmp(KbName(keyNamethr),'9') || strcmp(KbName(keyNamethr),'9(')==1
%                         sca
%                         return
%                     elseif strcmp(KbName(keyNamethr),'8') || strcmp(KbName(keyNamethr),'8*')==1
%                         break
%                     end
%                 elseif strcmp(KbName(keyNamePrac), '5')==1 || strcmp(KbName(keyNamethr),'5%')==1
%                     DrawFormattedText(screeninfo.win, textinfo.repinstr, 'center', 'center',screeninfo.fontcolour);
%                     Screen('DrawTexture', screeninfo.win, exampletargetr, [], gaborRectright);
%                     Screen('DrawTexture', screeninfo.win, exampletargetl, [], gaborRectleft);
%                     DrawFormattedText(screeninfo.win, orientation1, textleft(1),textleft(2),screeninfo.fontcolour);
%                     DrawFormattedText(screeninfo.win, orientation2, textright(1),textright(2),screeninfo.fontcolour);
%                     Screen('Flip', screeninfo.win);
%                     answered=1;
%                     KbStrokeWait;
%                 elseif strcmp(KbName(keyNamePrac),'2')==1 || strcmp(KbName(keyNamethr),'2@')==1
%                     startexp=1;
%                     answered=1;
%                 else
%                 end
%             end
%         end

        % Start
        i=1;

        counter.rhy=1;
        counter.int=1;
        counter.irr=1;


        for b=mixinfo.startblock:mixinfo.nBlockstot
            currentcondition=conditionmatrix(b);
            % Determine Current Condition
            if conditionmatrix(b)==1 % Rhythm
                blocklength = mixinfo.regularblength;
                textinfo.blocklabel= 'Rhythm';
                afcpass=0;
                afccount=0;
                if b<=3 %for the first rhythm block, initiate 2afc task
                    DrawFormattedText(screeninfo.win, textinfo.rhythmlong, 'center', 'center',screeninfo.fontcolour);
                    Screen('Flip', screeninfo.win);
                    KbStrokeWait;
                    Screen('DrawTexture', screeninfo.win, rhyimgtex);
                    if session_n==1
                    DrawFormattedText(screeninfo.win, 'Press any button to continue.', 'center', rhytextpos,screeninfo.fontcolour);
                    else
                    DrawFormattedText(screeninfo.win, 'Press any button to start the block.', 'center', rhytextpos,screeninfo.fontcolour);
                    end
                    Screen('Flip', screeninfo.win);
                    KbStrokeWait;
                    if twoafc && session_n==1
                        DrawFormattedText(screeninfo.win, textinfo.twoAFCinstr, 'center', 'center',screeninfo.fontcolour);
                        Screen('Flip', screeninfo.win);
                        KbStrokeWait;
                        while ~afcpass % Run 2AFC Task
                            afccount=afccount+1;
                            afceval(afccount)=twoafcfunction(currentcondition, screeninfo, visStiminfo, twoafcinfo, timeinfo);
                            if afceval(afccount)==1
                                DrawFormattedText(screeninfo.win, 'Correct', 'center', 'center',screeninfo.fontcolour);
                                Screen('Flip', screeninfo.win);
                                pause(0.5)
                            else
                                DrawFormattedText(screeninfo.win, 'Incorrect', 'center', 'center',screeninfo.fontcolour);
                                Screen('Flip', screeninfo.win);
                                pause(0.5)
                            end
                            if afccount>=twoafcinfo.minafceval % Start checking performance after the minimum of 2AFC trials
                                if sum(afceval(1,end-twoafcinfo.afcconseqcorr+1:end))==twoafcinfo.afcconseqcorr
                                    afcpass=1;
                                end
                            elseif afccount==twoafcinfo.maxafcrep % If Maxmimum trials are reached without reaching afcconseqcorr.
                                DrawFormattedText(screeninfo.win, 'Attention. Exceeded maximum trials. Repeat? (y/n)', 'center', 'center',screeninfo.fontcolour);
                                Screen('Flip', screeninfo.win);
                                responded=0;
                                while ~responded
                                    [keyDown, ~, keyName] = KbCheck;
                                    if keyDown
                                        if strcmp(KbName(keyName), 'y')==1
                                            responded=1;
                                        elseif strcmp(KbName(keyName), 'n')==1
                                            responded=1;
                                            afcpass=1;
                                        elseif strcmp(KbName(keyName),'e')==1
                                            DrawFormattedText(screeninfo.win,'Continue (c)? Exit (e)?', 'center', 'center', screeninfo.fontcolour);
                                            Screen('Flip', screeninfo.win);
                                            [~, keyNamethr, ~]=KbStrokeWait;
                                            if strcmp(KbName(keyNamethr),'e')
                                                sca
                                                return
                                            elseif strcmp(KbName(keyNamethr),'c')
                                                DrawFormattedText(screeninfo.win,'Ready? Press any button.', 'center', 'center', screeninfo.fontcolour);
                                                Screen('Flip', screeninfo.win);
                                                KbStrokeWait;
                                                break
                                            end
                                        end
                                    end
                                end
                            end
                        end
                        DrawFormattedText(screeninfo.win, ['Great. Let´s start with the real block. Please report the orientation of the target. \n\n' ...
                        'Use the timing to help you prepare for the right moment. The target will always be in line with the structure.'], 'center', 'center',screeninfo.fontcolour);
                        Screen('Flip', screeninfo.win);
                        KbStrokeWait;
                    end
                else % do not run 2afc again after the first block
                    DrawFormattedText(screeninfo.win, textinfo.rhythmshort, 'center', 'center',screeninfo.fontcolour);
                    Screen('Flip', screeninfo.win);
                    KbStrokeWait;
                end
                %             Screen('DrawTexture', screeninfo.win, rhyimgtex, [], illustrationRect);
                alltrials=alltrials1;
            elseif conditionmatrix(b)==2 % Interval
                currentcondition=conditionmatrix(b);
                blocklength = mixinfo.regularblength;
                textinfo.blocklabel = 'Interval';
                afcpass=0;
                afccount=0;
                if b<=3 %for the first rhythm block, initiate 2afc task
                    DrawFormattedText(screeninfo.win, textinfo.intervallong, 'center', 'center',screeninfo.fontcolour);
                    Screen('Flip', screeninfo.win);
                    KbStrokeWait;
                    Screen('DrawTexture', screeninfo.win, intimgtex);
                    if session_n==1
                    DrawFormattedText(screeninfo.win, 'Press any button to continue.', 'center', rhytextpos,screeninfo.fontcolour);
                    else
                    DrawFormattedText(screeninfo.win, 'Press any button to start the block.', 'center', rhytextpos,screeninfo.fontcolour);
                    end
                    Screen('Flip', screeninfo.win);
                    KbStrokeWait;
                    if twoafc
                        DrawFormattedText(screeninfo.win, textinfo.twoAFCinstr, 'center', 'center',screeninfo.fontcolour);
                        Screen('Flip', screeninfo.win);
                        KbStrokeWait;
                        while ~afcpass % Run 2AFC Task
                            afccount=afccount+1;
                            afceval(afccount)=twoafcfunction(currentcondition, screeninfo, visStiminfo, twoafcinfo, timeinfo);
                            if afceval(afccount)==1
                                DrawFormattedText(screeninfo.win, 'Correct', 'center', 'center',screeninfo.fontcolour);
                                Screen('Flip', screeninfo.win);
                                pause(0.5)
                            else
                                DrawFormattedText(screeninfo.win, 'Incorrect', 'center', 'center',screeninfo.fontcolour);
                                Screen('Flip', screeninfo.win);
                                pause(0.5)
                            end

                            if afccount>=twoafcinfo.minafceval % Start checking performance after the minimum of 2AFC trials
                                if sum(afceval(1,end-twoafcinfo.afcconseqcorr+1:end))==twoafcinfo.afcconseqcorr
                                    afcpass=1;
                                end
                            elseif afccount==twoafcinfo.maxafcrep % If Maxmimum trials are reached without reaching afcconseqcorr.
                                DrawFormattedText(screeninfo.win, 'Attention. Exceeded maximum trials. Repeat? (y/n)', 'center', 'center',screeninfo.fontcolour);
                                Screen('Flip', screeninfo.win);
                                responded=0;
                                while ~responded
                                    [keyDown, ~, keyName] = KbCheck;
                                    if keyDown
                                        if strcmp(KbName(keyName), 'y')==1
                                            responded=1;
                                        elseif strcmp(KbName(keyName), 'n')==1
                                            responded=1;
                                            afcpass=1;
                                        elseif strcmp(KbName(keyName),'e')==1
                                            DrawFormattedText(screeninfo.win,'Continue (c)? Exit (e)?', 'center', 'center', screeninfo.fontcolour);
                                            Screen('Flip', screeninfo.win);
                                            [~, keyNamethr, ~]=KbStrokeWait;
                                            if strcmp(KbName(keyNamethr),'e')
                                                sca
                                                return
                                            elseif strcmp(KbName(keyNamethr),'c')
                                                DrawFormattedText(screeninfo.win,'Ready? Press any button.', 'center', 'center', screeninfo.fontcolour);
                                                Screen('Flip', screeninfo.win);
                                                KbStrokeWait;
                                                break
                                            end
                                        end
                                    end
                                end
                            end
                        end
                         DrawFormattedText(screeninfo.win, ['Great. Let´s start with the real block. Please report the orientation of the target. \n\n' ...
                        'Use the timing to help you prepare for the right moment. The target will always be in line with the preceding time structure.'], 'center', 'center',screeninfo.fontcolour);
                         Screen('Flip', screeninfo.win);
                         KbStrokeWait;
                    end
                else
                    DrawFormattedText(screeninfo.win, textinfo.intervalshort, 'center', 'center',screeninfo.fontcolour);
                    Screen('Flip', screeninfo.win);
                    KbStrokeWait;
                end
                %             Screen('DrawTexture', screeninfo.win, intimgtex, [], illustrationRect);
                alltrials=alltrials2;
            elseif conditionmatrix(b)==3 % Irregular
                currentcondition=conditionmatrix(b);
                blocklength = mixinfo.irregularblength;
                textinfo.blocklabel = 'Irregular';
                alltrials=alltrials3;
                if b~=1
                DrawFormattedText(screeninfo.win, textinfo.irregularshort, 'center', 'center',screeninfo.fontcolour);
                Screen('Flip', screeninfo.win);
                KbStrokeWait;
                end
            end

            % Quick practice
            qpracticetrials(1,1:(nqp/2))=1;
            qpracticetrials(1,(nqp/2+1):nqp)=0;
            [m,n] = size(qpracticetrials);
            idx = randperm(n);
            y = qpracticetrials;
            y(1,idx) = qpracticetrials(1,:);
            qpracticetrials=y';
            qpcontrasts(1:nqp)=PracticeContrast;
            qpalltrials=[qpcontrasts' qpracticetrials];
            for qpidx=1:nqp
                [~, ~, OrRespEval, ~,~,~,~,~]=trialfunctionpro1(currentcondition,qpidx,0,screeninfo, visStiminfo, timeinfo, textinfo, qpalltrials,counter,1,0);
                if OrRespEval==1
                    DrawFormattedText(screeninfo.win, 'Correct', 'center', 'center',screeninfo.fontcolour);
                    Screen('Flip', screeninfo.win);
                    pause(1)
                else
                    DrawFormattedText(screeninfo.win, 'Incorrect', 'center', 'center',screeninfo.fontcolour);
                    Screen('Flip', screeninfo.win);
                    pause(1)
                end
            end
            DrawFormattedText(screeninfo.win, 'Practice Completed. Press any key to start the block.', 'center', 'center',screeninfo.fontcolour);
            Screen('Flip', screeninfo.win);
            KbStrokeWait;
            DrawFormattedText(screeninfo.win, 'Loading.', 'center', 'center',screeninfo.fontcolour);
            Screen('Flip', screeninfo.win);
            %Trials

            io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, 254); % start recording
            pause(5)
            io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, 0);

            for t=1:blocklength
                Screen('Flip', screeninfo.win,[],1);
                [RT, OrResp, OrRespEval, VisResp, ObjWarn, SubjWarn, FrameRate,IrrCueTimes]=trialfunctionpro1(currentcondition,t,0,screeninfo, visStiminfo, timeinfo, textinfo, alltrials,counter,0,0);

                % Save Trial Results
                resulttable(i,:)=table(subj, currentcondition, visStiminfo.maskintensity, b, t, FrameRate, RT, OrResp, OrRespEval, VisResp, ObjWarn, SubjWarn,visStiminfo.gaborpercent, 'VariableNames',{'Subject', ...
                    'Condition','Mask Intensity','Block','Trial','Frame Rate','Reaction Time', 'Orientation Reseponse', 'Correct/Incorrect', 'Visibility Response','Late Warning Obj','Late Warning Subj', ...
                    'Gabor Strength'});
                IrrCues(i,:)=IrrCueTimes;
                subresults.data=resulttable;
                subresults.irrcuetimes=IrrCues;
                save([fileNameAll '.mat'], 'subresults')
                i=i+1;

                % Update counters: (only here in case one trial gets interrupted and
                % repeated)
                if currentcondition==1
                    counter.rhy=counter.rhy+1;
                elseif currentcondition==2
                    counter.int=counter.int+1;
                elseif currentcondition==3
                    counter.irr=counter.irr+1;
                end

                %             % give chance to interrupt speedrun
                %             if screeninfo.speedrun
                %                 for i=1:120
                %                     [keyIsDown, ~, keyCode, ~] = KbCheck;
                %                     if keyIsDown
                %                         if strcmp(KbName(keyCode),'e')==1
                %                             DrawFormattedText(screeninfo.win,'Continue (c)? Exit (e)?', 'center', 'center', screeninfo.fontcolour);
                %                             Screen('Flip', screeninfo.win);
                %                             [~, keyNamethr, ~]=KbStrokeWait;
                %                             if strcmp(KbName(keyNamethr),'e')
                %                                 sca
                %                                 return
                %                             elseif strcmp(KbName(keyNamethr),'c')
                %                                 continue
                %                             end
                %                         else
                %                         end
                %                     else
                %                     end
                %                 end
                %             else
                %             end
                %       KbStrokeWait;
            end


            % Calculate block performance
            blockperf=resulttable{find(resulttable{:,'Block'}==b),'Correct/Incorrect'}; % performance results for latest block
            blockaccuracy=mean(blockperf,'omitnan'); % mean accuracy
            totalperf=resulttable{:,'Correct/Incorrect'};
            totalaccuracy= mean(totalperf,'omitnan'); % total accuracy of all trials SO FAR
            accuracy=[blockaccuracy totalaccuracy];
            
            % DrawFormattedText(screeninfo.win, 'Loading.', 'center', 'center',screeninfo.fontcolour);
            % Screen('Flip', screeninfo.win);
            pause(5)
            io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, 255); % Stop recording
            pause(1)
            io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, 0);

            % Pause Between Blocks (and disaplay accuracy)
            if b<3
                breaktext=sprintf('Block %i/%i completed. Please take a break.\n \n\n \n Wait for the experimenter to start the next block.', b,(mixinfo.nBlockstot));
            elseif b==mixinfo.nBlockstot % after last block
                breaktext=sprintf('Block %i/%i completed. This is the end of the experiment.\n \n\n \n Please wait for the experimenter.', b,(mixinfo.nBlockstot));
            else
                breaktext=sprintf('Block %i/%i completed. Please take a break.\n \n\n \n Press any button to start the next block.', b,(mixinfo.nBlockstot));
            end
            DrawFormattedText(screeninfo.win, breaktext, 'center', 'center',screeninfo.fontcolour);
            performancetext=sprintf('B%.2f24753 /n T%.2f45264', accuracy(1),accuracy(2)); % block and total
            [~,~,heightacc,lengthacc]=Screen('TextBounds', screeninfo.win, performancetext); % get bounds
            DrawFormattedText(screeninfo.win, performancetext, screeninfo.axisx(1), screeninfo.axisy(2)-heightacc*1.2,screeninfo.fontcolour);
            Screen('Flip', screeninfo.win);

            % give option to adjust difficulty after first block
            if b==1
                finishadjusting2=0;
                adjustedafterblock1=0;
                adjusttext=sprintf('This is the end of the first block.\n\n Please wait for the experimenter \n\n \n\n%.2f \n Adjust? (y-2/n-3)', accuracy(1));
                DrawFormattedText(screeninfo.win,adjusttext, 'center', 'center',screeninfo.fontcolour);
                Screen('Flip', screeninfo.win);
                while finishadjusting2==0
                    [~, keyNamethr, ~]=KbStrokeWait;
                    if strcmp(KbName(keyNamethr),'9')==1 || strcmp(KbName(keyNamethr),'9(')==1
                        DrawFormattedText(screeninfo.win,'Continue (8)? Exit (9)?', 'center', 'center', screeninfo.fontcolour);
                        Screen('Flip', screeninfo.win);
                        [~, keyNamethr, ~]=KbStrokeWait;
                        if strcmp(KbName(keyNamethr),'9') || strcmp(KbName(keyNamethr),'9(')==1
                            finishadjusting2=1;
                            sca
                            return
                        elseif strcmp(KbName(keyNamethr),'8') || strcmp(KbName(keyNamethr),'8*')==1
                            continue
                        end
                    elseif strcmp(KbName(keyNamethr), '2')==1 || strcmp(KbName(keyNamethr),'2@')==1
                        Screen('Flip', screeninfo.win);
                        while ~finishadjusting2
                            contrasttext='Current intensity: %.2f \n\n Type new intensity (0.01-0.99):';
                            newintensity=adjustintensity(screeninfo, visStiminfo.gaborpercent);
                            ThresholdReachedText=sprintf('New intensity: %.2f \n\n Correct? (y-2/n-3)',newintensity);
                            DrawFormattedText(screeninfo.win,ThresholdReachedText, 'center', 'center', screeninfo.fontcolour);
                            Screen('Flip', screeninfo.win);
                            [~, keyNamethr, ~]=KbStrokeWait;
                            if strcmp(KbName(keyNamethr),'9')==1 || strcmp(KbName(keyNamethr),'9(')==1
                                DrawFormattedText(screeninfo.win,'Continue (8)? Exit (9)?', 'center', 'center', screeninfo.fontcolour);
                                Screen('Flip', screeninfo.win);
                                [~, keyNamethr, ~]=KbStrokeWait;
                                if strcmp(KbName(keyNamethr),'9') || strcmp(KbName(keyNamethr),'9(')==1
                                    sca
                                    return
                                elseif strcmp(KbName(keyNamethr),'8') || strcmp(KbName(keyNamethr),'8*')==1
                                    continue 
                                end
                            elseif strcmp(KbName(keyNamethr), '2')==1 || strcmp(KbName(keyNamethr),'2@')==1 %continue
                                visStiminfo.gaborpercent=newintensity;
                                finishadjusting2=1;
                                adjustedafterblock1=1;
                            elseif strcmp(KbName(keyNamethr),'3')==1 || strcmp(KbName(keyNamethr),'3#')==1 % type threshold again
                            end
                        end
                    elseif strcmp(KbName(keyNamethr),'3')==1 || strcmp(KbName(keyNamethr),'3#')==1
                        finishadjusting2=1;
                    else
                    end
                end

            end



            % Display on second monitor
            %             figuretitle=sprintf('Accuracy in Block %i: %.2f', b, blockaccuracy);
            %             figure;
            %             plot(accuracy,'o','LineWidth',2)
            %             title(figuretitle)
            %             ylabel('Accuracy in %')
            %             xlim([0 3])
            %             ylim([-0.1 1.1])
            %             xticks(1:2)
            %             xticklabels({'Block Accuracy','Total Mean Accuracy'})


            if screeninfo.speedrun
                [~, keyNamethr, ~]=KbStrokeWait;
                if strcmp(KbName(keyNamethr),'e')==1
                    DrawFormattedText(screeninfo.win,'Continue (c)? Exit (e)?', 'center', 'center', screeninfo.fontcolour);
                    Screen('Flip', screeninfo.win);
                    [~, keyNamethr, ~]=KbStrokeWait;
                    if strcmp(KbName(keyNamethr),'e')
                        sca
                        return
                    elseif strcmp(KbName(keyNamethr),'c')
                        continue
                    end
                else
                end
            end
        end
       
        % Save
        subresults.data=resulttable;
        subresults.conditions=conditionmatrix;
        save([fileNameAll '.mat'], 'subresults')

        % End Screen
        DrawFormattedText(screeninfo.win, textinfo.endtext, 'center', 'center', screeninfo.fontcolour);
        Screen('Flip', screeninfo.win);
        KbStrokeWait;

        % Clear screen
        sca;
        ShowCursor;

    elseif screeninfo.scopetest % If scope test, run individual trials of different conditions, as chosen in the initial menu
        % add scope test here
        finishscope=0;
        while ~finishscope
            condchosen=0;
            while ~condchosen
                DrawFormattedText(screeninfo.win, 'Oscilloscope Test. Choose condition. \n\n Rhythm(1) \nInterval (2) \nIrregular (3)\n Exit (e)', 'center', 'center', screeninfo.fontcolour);
                Screen('Flip', screeninfo.win);

                [~, keyNamethr, ~]=KbStrokeWait;
                if strcmp(KbName(keyNamethr),'e')
                    finishscope=1;
                    sca
                    return
                elseif strcmp(KbName(keyNamethr),'1')
                    currentcondition=1;
                    condchosen=1;
                elseif strcmp(KbName(keyNamethr),'2')
                    currentcondition=2;
                    condchosen=1;
                elseif strcmp(KbName(keyNamethr),'3')
                    currentcondition=3;
                    condchosen=1;
                end
            end
            % Run trial
            [~, ~, ~, ~,~,~,~,~]=trialfunctionpro1(currentcondition,1,0,screeninfo, visStiminfo, timeinfo, textinfo, 0, 0,0,0);
        end
    end
catch
    sca;
    ShowCursor;
    psychrethrow(psychlasterror);
end

%% Trial Matrix Function
function [alltrialsmixed]=fulltrialmatrix(mixinfo,irregular)

cont=1:10;
orientations=0:1;
irrtargetbins=1:5;


if irregular
    % create a batch of 10 trials with 2x5 time bins
    binsOneUnit=[irrtargetbins irrtargetbins]';

    % repeat for each batch of 10, but with shifted time bins to counterbalance
    binsAll=[];
    for i=0:4
        binsAll=[binsAll; mod(binsOneUnit+i,5)];
    end
    binsAll(binsAll==0)=5; % finish shifting

    % add contrast levelsall
    designMat=[kron(ones(1,5), cont)' binsAll];

    % randomizing the batches of 10 within block (diff for each participant)
    order=randperm(5);

    % Add orientation and shuffle
    orient=kron(orientations, ones(1,5))';
    orient=orient(randperm(10));  % forgot why we would need to shuffle here. remove?

    % counterbalance orientation for each second batch of 10
    block1=[];
    for j=1:5

        take=designMat((order(j)-1)*10+1:order(j)*10,:);
        if mod(j,2)==0
            orientCorrect=orient;
        else
            orientCorrect=1-orient;
        end
        block1=[block1; [take orientCorrect]];

    end

    % Create second block and concatenate for all irregular trials
    block2=block1;
    block2(:,3)=1-block2(:,3);
    fullDesigntemp=[block1; block2; block1; block2];

    % switch columns to match regular condition matrices
    fullDesign=fullDesigntemp;
    fullDesign(:,2)=fullDesigntemp(:,3);
    fullDesign(:,3)=fullDesigntemp(:,2);
    
    %
    catchtrials=zeros(mixinfo.nCatchSampleirr,2); % add any value in third row for catch trials
    catchtrials=[catchtrials ones(mixinfo.nCatchSampleirr,1)];
else
    batch1=[kron([1 1], cont)' kron(ones(10,1), orientations')];
    batch1(11:20,2)=1-batch1(11:20,2);
    fullDesign=[];
    catchtrials=zeros(mixinfo.nCatchSample,2);
    for i=1:6 %hard coded, should be parameterized (needs to lead to total amount of trials per condition
        fullDesign=[fullDesign; batch1];
    end
end


%% Shuffle trials within each batch of 10 and according ot constraints

nbatches=length(fullDesign)/length(cont); % how many batches of 10 in total condition?
alltrialsmixed=[];

contrastlevels=length(cont); % initialize counter
start=1;
prevshuffledbatch=[];

for nb=1:nbatches
    % Select current batch
    currbatch=fullDesign(start:contrastlevels,:);

 currbatch=[currbatch; catchtrials]; % combine target and catch trials

    % Shuffle
    conditionsMet=0;
    irrmet=0;
    batchlength=length(currbatch);

    while ~conditionsMet
        shuffledBatch=currbatch(randperm(batchlength), :); % shuffle

        % create an extended 'shuffled batch' including the last two trials
        % of the previous block so all following rules also cover the
        % merging point of different batches
        if nb>1
            shuffledBatchext=[prevshuffledbatch((end-1:end),:);shuffledBatch];
        else
            shuffledBatchext=shuffledBatch;
        end

        % Rules
        catchAtStart=shuffledBatch(1,1)==0; % first trial not a catch trial
        catchConsec=max(diff([0; find(shuffledBatchext(:,1)~=0); (size(shuffledBatchext,1)+1)])-1); % maximum catch trialltrialsals in a row
        leftConsec=max(diff([0; find(shuffledBatchext(:,2)~=0); (size(shuffledBatchext,2)+1)])-1); % orientation consequtive
        rightConsec=max(diff([0; find(shuffledBatchext(:,2)~=1); (size(shuffledBatchext,2)+1)])-1); % orientation consequtive

        % Extra rules for irregular condition
        if irregular
            shortConsec=max(diff([0; find(shuffledBatchext(:,3)>=3); (size(shuffledBatchext,2)+1)])-1); % short time bin consequtive
            longConsec=max(diff([0; find(shuffledBatchext(:,3)<=3); (size(shuffledBatchext,2)+1)])-1); % short time bin consequtive
        else
        end

        % All rules fulfilled? Then finish, otherwise repeat
        if irregular
            if shortConsec <= mixinfo.maxshortrep && longConsec <= mixinfo.maxlongrep
                irrmet=1;
            else
                irrmet=0;
            end
        else
            irrmet=1;
        end

        if catchAtStart==0 && catchConsec<=mixinfo.maxconscatch && leftConsec<=mixinfo.maxorientrep && rightConsec<=mixinfo.maxorientrep...
                && irrmet==1
            conditionsMet=1;
            alltrialsmixed=[alltrialsmixed;shuffledBatch];
        else
            conditionsMet=0;
        end
    end
    start=start+mixinfo.nContrasts;
    contrastlevels=contrastlevels+mixinfo.nContrasts;
    prevshuffledbatch=shuffledBatch;
end
end

%% Threshold Function
function [OrRespEval, VisResp,repeat]=thresholdfunctionpro1(currcontrast,screeninfo, visStiminfo, timeinfo, textinfo)

% Choose orientation
thrtrialorient = randi([0,1],1);
if thrtrialorient == 0
    visStiminfo.gaborcurr=visStiminfo.gaborRight; %Right
else
    visStiminfo.gaborcurr=visStiminfo.gaborLeft; %left
end

% Create Gabor embedded in mask
for i=1:timeinfo.TarDurFr
    if ~screeninfo.scopetest
        GaborAdjusted=currcontrast*visStiminfo.gaborcurr;
        randMask=(rand(visStiminfo.rectSize)-0.5)*2;
        targetStim=0.5+0.5*(visStiminfo.maskintensity*randMask+GaborAdjusted*visStiminfo.gaborpercent); % input: between -1 and 1, output between 0 and 1)
        visStiminfo.currTarget(i)=Screen('MakeTexture', screeninfo.win, targetStim);
    else
        noiseimg=0.01+visStiminfo.maskintensity*(rand(visStiminfo.rectSize)-0.5);
        % make cue area of noise tex white
        xfrom=(length(visStiminfo.noiseimg)-visStiminfo.baseRect(3))/2;
        xto=xfrom+visStiminfo.baseRect(3);
        noiseimg(xfrom:xto,xfrom:xto)=visStiminfo.rectColorCue;
        visStiminfo.currTarget(i)=Screen('MakeTexture', screeninfo.win, noiseimg);
    end
end

% Shuffle noise and initialize counter
visStiminfo.noisetex=visStiminfo.noisetex(randperm(length(visStiminfo.noisetex)));
cnoise=1;

% Choose jittered pre target time
pretarjitt = randi(3);

% Present Noise
for fixtime=1:timeinfo.thrFr(pretarjitt)
    Screen('DrawTexture', screeninfo.win, visStiminfo.noisetex(cnoise), [], [], 0);
    Screen('Flip', screeninfo.win);
    cnoise=cnoise+1;
end

% Present Target
for TarTime=1:timeinfo.TarDurFr
    cnoise=cnoise+1;
    Screen('DrawTexture', screeninfo.win, visStiminfo.currTarget(TarTime));
    Screen('Flip', screeninfo.win);
end

for QuestInt=1:timeinfo.TargQFr
    cnoise=cnoise+1;
    Screen('DrawTexture', screeninfo.win, visStiminfo.noisetex(cnoise), [], [], 0);
    Screen('Flip', screeninfo.win);
end

% Screen('DrawLines', screeninfo.win, visStiminfo.fixCoords,...
%     visStiminfo.lineWidthPix, screeninfo.black, [screeninfo.xCenter screeninfo.yCenter], 2);
% Screen('Flip', screeninfo.win);

% Request Orientation Reponse
DrawFormattedText(screeninfo.win, textinfo.OrientRespPrompt, 'center', 'center', screeninfo.fontcolour);
Screen('Flip', screeninfo.win);
orresponded=0;
visresponded=0;
repeat=0;
while ~orresponded
    [keyDown, ~, keyName] = KbCheck;
    if keyDown
        if strcmp(KbName(keyName), textinfo.KeyLeft)==1
            orresponded=1;
            if thrtrialorient==1 % left
                OrRespEval=1;
            else
                OrRespEval=0;
            end
        elseif strcmp(KbName(keyName), textinfo.KeyRight)==1
            orresponded=1;
            if thrtrialorient==0 %Right
                OrRespEval=1;
            else
                OrRespEval=0;
            end
        elseif strcmp(KbName(keyName),'e')==1
            DrawFormattedText(screeninfo.win,'Continue (c)? Exit (e)?', 'center', 'center', screeninfo.fontcolour);
            Screen('Flip', screeninfo.win);
            [~, keyNamethr, ~]=KbStrokeWait;
            if strcmp(KbName(keyNamethr),'e')
                sca
                return
            elseif strcmp(KbName(keyNamethr),'c')
                OrRespEval=NaN;
                VisResp=NaN; % do not stop threshold but instead repeat the interrupted round
                visresponded=1;
                repeat=1;
                DrawFormattedText(screeninfo.win,'Ready? Press any button.', 'center', 'center', screeninfo.fontcolour);
                Screen('Flip', screeninfo.win);
                KbStrokeWait
                break
            end
        end
    end
end

% Request Visibility Response
Screen('Flip', screeninfo.win);
pause(0.2)
while ~visresponded
    DrawFormattedText(screeninfo.win, textinfo.VisRespPrompt, 'center', 'center', screeninfo.fontcolour);
    Screen('Flip', screeninfo.win);
    [~, keyNamethr, ~]=KbStrokeWait;
    if strcmp(KbName(keyNamethr), textinfo.PAS1)==1
        visresponded=1;
        VisResp=0;
    elseif strcmp(KbName(keyNamethr), textinfo.PAS2)==1
        visresponded=1;
        VisResp=1;
    elseif strcmp(KbName(keyNamethr), textinfo.PAS3)==1
        visresponded=1;
        VisResp=2;
    elseif strcmp(KbName(keyNamethr), textinfo.PAS4)==1
        visresponded=1;
        VisResp=3;
    elseif strcmp(KbName(keyNamethr),'e')==1
        DrawFormattedText(screeninfo.win,'Continue (c)? Exit (e)?', 'center', 'center', screeninfo.fontcolour);
        Screen('Flip', screeninfo.win);
        [~, keyNamethr, ~]=KbStrokeWait;
        if strcmp(KbName(keyNamethr),'e')
            sca
            return
        elseif strcmp(KbName(keyNamethr),'c')
            VisResp=NaN; % don't record results for interrupted trials
            OrRespEval=NaN;
            DrawFormattedText(screeninfo.win,'Ready? Press any button.', 'center', 'center', screeninfo.fontcolour);
            Screen('Flip', screeninfo.win);
            KbStrokeWait
            break
        end
    else
    end
end
Screen('Flip', screeninfo.win);
end

%% Function: 2AFC Timing Task
function [respeval]=twoafcfunction(currentcondition, screeninfo, visStiminfo, twoafcinfo, timeinfo)

% Shuffle textures and initialize counters
cnoise=1;
ccue=1;
cwarn=1;
visStiminfo.noisetex=visStiminfo.noisetex(randperm(length(visStiminfo.noisetex)));
visStiminfo.warntex=visStiminfo.warntex(randperm(length(visStiminfo.warntex)));
visStiminfo.cuetex=visStiminfo.cuetex(randperm(length(visStiminfo.cuetex)));

% Timing:
if currentcondition==2 % for interval choose a jittered between-interval time
    bintTime = randi([1400,1800])/1000; %%% use parameter %% % currently around 1600 to match total rhythm length
    timeinfo.bintTimefr=round(bintTime/screeninfo.ifi);
end

% Create vertical Gabor embedded in mask for 2afc
for i=1:(timeinfo.TarDurFr*3)
        randMask=(rand(visStiminfo.rectSize)-0.5)*2;
        afctargetStim=0.5+0.5*(visStiminfo.maskintensity*randMask+visStiminfo.gabor2afc*1); % input: between -1 and 1, output between 0 and 1)
        visStiminfo.afcTarget(i)=Screen('MakeTexture', screeninfo.win, afctargetStim);
end


% Choose randomly which interval is correct
currcorrint=twoafcinfo.corrint(randi(2));

% Assign timing accordingly:
if currcorrint==1
    targettiming1=twoafcinfo.correct_timing_fr;
    targettiming2=twoafcinfo.offbeat_timing_fr(randi(2)); %randomly choose if it is too short or too long
elseif currcorrint==2
    targettiming1=twoafcinfo.offbeat_timing_fr(randi(2));
    targettiming2=twoafcinfo.correct_timing_fr;
end


% Present First Interval
intervallabel='1';
if currentcondition==1 % If rhythm
    % Present Rhythm

    for fixtime=1:timeinfo.ITIFr(randi(length(timeinfo.ITIFr)))
        Screen('DrawTexture', screeninfo.win, visStiminfo.noisetex(cnoise), [], [], 0);
        DrawFormattedText(screeninfo.win,'1','center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
        Screen('DrawLines', screeninfo.win, visStiminfo.fixCoords,...
            visStiminfo.lineWidthPix, screeninfo.black, [screeninfo.xCenter screeninfo.yCenter], 2);
        Screen('Flip', screeninfo.win);
        cnoise=cnoise+1;
    end

    for rhythms=1:3
        for CueDura=1:timeinfo.CueDurFr
            DrawFormattedText(screeninfo.win,'1','center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
            Screen('DrawTexture', screeninfo.win, visStiminfo.cuetex(ccue), [], [], 0);
            Screen('Flip', screeninfo.win);
            ccue=ccue+1;
        end

        for rhythISI=1:twoafcinfo.correct_timing_fr
            DrawFormattedText(screeninfo.win,'1','center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
            Screen('DrawTexture', screeninfo.win, visStiminfo.noisetex(cnoise), [], [], 0);
            Screen('Flip', screeninfo.win);
            cnoise=cnoise+1;
        end
    end

    % Present Target Interval
    for WarnTime=1:timeinfo.CueDurFr
        DrawFormattedText(screeninfo.win,'1','center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
        Screen('DrawTexture', screeninfo.win, visStiminfo.warntex(cwarn), [], [], 0);
        Screen('Flip', screeninfo.win);
        cwarn=cwarn+1;
    end

    for rhythISI=1:targettiming1
        DrawFormattedText(screeninfo.win,'1','center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
        Screen('DrawTexture', screeninfo.win, visStiminfo.noisetex(cnoise), [], [], 0);
        Screen('Flip', screeninfo.win);
        cnoise=cnoise+1;
    end

elseif currentcondition==2
    intervallabel='1';
    for fixtime=1:timeinfo.ITIFr(randi(length(timeinfo.ITIFr)))
        Screen('DrawTexture', screeninfo.win, visStiminfo.noisetex(cnoise), [], [], 0);
        DrawFormattedText(screeninfo.win,intervallabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
        Screen('DrawLines', screeninfo.win, visStiminfo.fixCoords,...
            visStiminfo.lineWidthPix, screeninfo.black, [screeninfo.xCenter screeninfo.yCenter], 2);
        Screen('Flip', screeninfo.win);
        cnoise=cnoise+1;
    end

    % for intervals=1:2
    for CueDura=1:timeinfo.CueDurFr
        DrawFormattedText(screeninfo.win,intervallabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
        Screen('DrawTexture', screeninfo.win, visStiminfo.cuetex(ccue), [], [], 0);
        Screen('Flip', screeninfo.win);
        ccue=ccue+1;
    end

    for intervalISI=1:twoafcinfo.correct_timing_fr
        DrawFormattedText(screeninfo.win,intervallabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
        Screen('DrawTexture', screeninfo.win, visStiminfo.noisetex(cnoise), [], [], 0);
        Screen('Flip', screeninfo.win);
        cnoise=cnoise+1;
    end

    for CueDura=1:timeinfo.CueDurFr
        DrawFormattedText(screeninfo.win,intervallabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
        Screen('DrawTexture', screeninfo.win, visStiminfo.cuetex(ccue), [], [], 0);
        Screen('Flip', screeninfo.win);
        ccue=ccue+1;
    end

    for intintISI=1:timeinfo.bintTimefr
        DrawFormattedText(screeninfo.win,intervallabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
        Screen('DrawTexture', screeninfo.win, visStiminfo.noisetex(cnoise), [], [], 0);
        Screen('Flip', screeninfo.win);
        cnoise=cnoise+1;
    end

    % Present Target Interval
    for WarnTime=1:timeinfo.CueDurFr
        DrawFormattedText(screeninfo.win,intervallabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
        Screen('DrawTexture', screeninfo.win, visStiminfo.warntex(cwarn), [], [], 0);
        Screen('Flip', screeninfo.win);
        cwarn=cwarn+1;
    end

    for intervalISI=1:targettiming1
        DrawFormattedText(screeninfo.win,intervallabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
        Screen('DrawTexture', screeninfo.win, visStiminfo.noisetex(cnoise), [], [], 0);
        Screen('Flip', screeninfo.win);
        cnoise=cnoise+1;
    end
end


% Present Target
for TarTime=1:(timeinfo.TarDurFr*3)
    DrawFormattedText(screeninfo.win,intervallabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
    Screen('DrawTexture', screeninfo.win, visStiminfo.afcTarget(TarTime));
    Screen('Flip', screeninfo.win);
    cwarn=cwarn+1;
end

for QuestInt=1:timeinfo.TargQFr
    DrawFormattedText(screeninfo.win,intervallabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
    Screen('DrawTexture', screeninfo.win, visStiminfo.noisetex(cnoise), [], [], 0);
    Screen('Flip', screeninfo.win);
    cnoise=cnoise+1;
end

DrawFormattedText(screeninfo.win,intervallabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
Screen('DrawLines', screeninfo.win, visStiminfo.fixCoords,...
    visStiminfo.lineWidthPix, screeninfo.black, [screeninfo.xCenter screeninfo.yCenter], 2);
Screen('Flip', screeninfo.win);

pause(1)


% Present Second Interval

% Shuffle textures and initialize counters
cnoise=1;
ccue=1;
cwarn=1;
visStiminfo.noisetex=visStiminfo.noisetex(randperm(length(visStiminfo.noisetex)));
visStiminfo.warntex=visStiminfo.warntex(randperm(length(visStiminfo.warntex)));
visStiminfo.cuetex=visStiminfo.cuetex(randperm(length(visStiminfo.cuetex)));
intervallabel='2';
if currentcondition==1 % If rhythm
    % Present Rhythm
    for fixtime=1:timeinfo.ITIFr(randi(length(timeinfo.ITIFr)))
        Screen('DrawTexture', screeninfo.win, visStiminfo.noisetex(cnoise), [], [], 0);
        DrawFormattedText(screeninfo.win,intervallabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
        Screen('DrawLines', screeninfo.win, visStiminfo.fixCoords,...
            visStiminfo.lineWidthPix, screeninfo.black, [screeninfo.xCenter screeninfo.yCenter], 2);
        Screen('Flip', screeninfo.win);
        cnoise=cnoise+1;
    end

    for rhythms=1:3
        for CueDura=1:timeinfo.CueDurFr
            DrawFormattedText(screeninfo.win,intervallabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
            Screen('DrawTexture', screeninfo.win, visStiminfo.cuetex(ccue), [], [], 0);
            Screen('Flip', screeninfo.win);
            ccue=ccue+1;
        end

        for rhythISI=1:twoafcinfo.correct_timing_fr
            DrawFormattedText(screeninfo.win,intervallabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
            Screen('DrawTexture', screeninfo.win, visStiminfo.noisetex(cnoise), [], [], 0);
            Screen('Flip', screeninfo.win);
            cnoise=cnoise+1;
        end
    end

    % Present Target Interval
    for WarnTime=1:timeinfo.CueDurFr
        DrawFormattedText(screeninfo.win,intervallabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
        Screen('DrawTexture', screeninfo.win, visStiminfo.warntex(cwarn), [], [], 0);
        Screen('Flip', screeninfo.win);
        cwarn=cwarn+1;
    end

    for rhythISI=1:targettiming2
        DrawFormattedText(screeninfo.win,intervallabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
        Screen('DrawTexture', screeninfo.win, visStiminfo.noisetex(cnoise), [], [], 0);
        Screen('Flip', screeninfo.win);
        cnoise=cnoise+1;
    end

elseif currentcondition==2
    for fixtime=1:timeinfo.ITIFr(randi(length(timeinfo.ITIFr)))
        Screen('DrawTexture', screeninfo.win, visStiminfo.noisetex(cnoise), [], [], 0);
        DrawFormattedText(screeninfo.win,intervallabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
        Screen('DrawLines', screeninfo.win, visStiminfo.fixCoords,...
            visStiminfo.lineWidthPix, screeninfo.black, [screeninfo.xCenter screeninfo.yCenter], 2);
        Screen('Flip', screeninfo.win);
        cnoise=cnoise+1;
    end

    % for intervals=1:2
    for CueDura=1:timeinfo.CueDurFr
        DrawFormattedText(screeninfo.win,intervallabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
        Screen('DrawTexture', screeninfo.win, visStiminfo.cuetex(ccue), [], [], 0);
        Screen('Flip', screeninfo.win);
        ccue=ccue+1;
    end

    for intervalISI=1:twoafcinfo.correct_timing_fr
        DrawFormattedText(screeninfo.win,intervallabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
        Screen('DrawTexture', screeninfo.win, visStiminfo.noisetex(cnoise), [], [], 0);
        Screen('Flip', screeninfo.win);
        cnoise=cnoise+1;
    end

    for CueDura=1:timeinfo.CueDurFr
        DrawFormattedText(screeninfo.win,intervallabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
        Screen('DrawTexture', screeninfo.win, visStiminfo.cuetex(ccue), [], [], 0);
        Screen('Flip', screeninfo.win);
        ccue=ccue+1;
    end

    for intintISI=1:timeinfo.bintTimefr
        DrawFormattedText(screeninfo.win,intervallabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
        Screen('DrawTexture', screeninfo.win, visStiminfo.noisetex(cnoise), [], [], 0);
        Screen('Flip', screeninfo.win);
        cnoise=cnoise+1;
    end

    % Present Target Interval
    for WarnTime=1:timeinfo.CueDurFr
        DrawFormattedText(screeninfo.win,intervallabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
        Screen('DrawTexture', screeninfo.win, visStiminfo.warntex(cwarn), [], [], 0);
        Screen('Flip', screeninfo.win);
        cwarn=cwarn+1;
    end

    for intervalISI=1:targettiming2
        DrawFormattedText(screeninfo.win,intervallabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
        Screen('DrawTexture', screeninfo.win, visStiminfo.noisetex(cnoise), [], [], 0);
        Screen('Flip', screeninfo.win);
        cnoise=cnoise+1;
    end
end


% Present Target
for TarTime=1:(timeinfo.TarDurFr*3)
    DrawFormattedText(screeninfo.win,intervallabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
    Screen('DrawTexture', screeninfo.win, visStiminfo.afcTarget(TarTime));
    Screen('Flip', screeninfo.win);
    cwarn=cwarn+1;
end

for QuestInt=1:timeinfo.TargQFr
    DrawFormattedText(screeninfo.win,intervallabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
    Screen('DrawTexture', screeninfo.win, visStiminfo.noisetex(cnoise), [], [], 0);
    Screen('Flip', screeninfo.win);
    cnoise=cnoise+1;
end

DrawFormattedText(screeninfo.win,intervallabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
Screen('DrawLines', screeninfo.win, visStiminfo.fixCoords,...
    visStiminfo.lineWidthPix, screeninfo.black, [screeninfo.xCenter screeninfo.yCenter], 2);
Screen('Flip', screeninfo.win);

pause(1)

% Request Response
DrawFormattedText(screeninfo.win, 'Which presentation matched the timing? \n\n 1 or 2?', 'center', 'center', screeninfo.fontcolour);
Screen('Flip', screeninfo.win);
responded=0;
while ~responded
    [keyDown, ~, keyName] = KbCheck;
    if keyDown
        if strcmp(KbName(keyName), '1')==1
            responded=1;
            if currcorrint==1 % left
                respeval=1;
            else
                respeval=0;
            end
        elseif strcmp(KbName(keyName), '2')==1
            responded=1;
            if currcorrint==2 %Right
                respeval=1;
            else
                respeval=0;
            end
        elseif strcmp(KbName(keyName),'e')==1
            DrawFormattedText(screeninfo.win,'Continue (c)? Exit (e)?', 'center', 'center', screeninfo.fontcolour);
            Screen('Flip', screeninfo.win);
            [~, keyNamethr, ~]=KbStrokeWait;
            if strcmp(KbName(keyNamethr),'e')
                sca
                return
            elseif strcmp(KbName(keyNamethr),'c')
                DrawFormattedText(screeninfo.win,'Ready? Press any button.', 'center', 'center', screeninfo.fontcolour);
                Screen('Flip', screeninfo.win);
                KbStrokeWait
                break
            end
        end
    end
end

end
%% Trial function
function [RT, OrResp, OrRespEval, VisResp,warningobj,subwarning,curfrate,irrISI]=trialfunctionpro1(currentcondition,t,firstprac,screeninfo, visStiminfo, timeinfo, textinfo, alltrials,counter,qp,justtar)

% Get current frame rate
curfrate=Screen('NominalFrameRate', screeninfo.win);

if timeinfo.recalcframes
    screeninfo.ifi = Screen('GetFlipInterval', screeninfo.win);

    % Recalculate frames with current timing info
    timeinfo.rhyISIfr=round(timeinfo.rhyISI/screeninfo.ifi);
    timeinfo.intISIfr=round(timeinfo.intISI/screeninfo.ifi);
    timeinfo.CueDurFr=round(timeinfo.CueDur/screeninfo.ifi);
    timeinfo.TarDurFr=round(timeinfo.TarDur/screeninfo.ifi);
    timeinfo.TargQFr=round(timeinfo.TargetQuestInt/screeninfo.ifi);
    timeinfo.ITIFr=round(timeinfo.ITI/screeninfo.ifi);
    timeinfo.thrFr=round(timeinfo.prethr/screeninfo.ifi);
    timeinfo.irrtartimefr=round(timeinfo.irrtartime/screeninfo.ifi);
end

% Choose counter
if currentcondition==1 && qp==0 && ~screeninfo.scopetest % don't increase counters for quick practice or scope tests
    count=counter.rhy;
elseif currentcondition==2 && qp==0 && ~screeninfo.scopetest
    count=counter.int;
elseif currentcondition==3 && qp==0 && ~screeninfo.scopetest
    count=counter.irr;
elseif qp
    count=t;
else
    count=1;
end

% Quick practice?

if currentcondition==0
    blocklabel='Practice';
elseif qp
    blocklabel=append(textinfo.blocklabel,' - ','Practice');
elseif ~currentcondition==0 && ~screeninfo.scopetest
    blocklabel=textinfo.blocklabel;
elseif screeninfo.scopetest
    blocklabel=' ';
end

% Choose orientation
 if ~screeninfo.scopetest
    trialorient=alltrials(count,2);
    if trialorient == 0
        visStiminfo.gaborcurr=visStiminfo.gaborRight; %Right
    else
        visStiminfo.gaborcurr=visStiminfo.gaborLeft; %left
    end
 else
     trialorient=randi(2);
 end

% Shuffle textures and initialize counters
cnoise=1;
ccue=1;
cwarn=1;
visStiminfo.warntex=visStiminfo.warntex(randperm(length(visStiminfo.warntex)));
visStiminfo.cuetex=visStiminfo.cuetex(randperm(length(visStiminfo.cuetex)));

% Choose interval between target and first question
preqtime=timeinfo.preqjitter(randi(3));
pretarjitt = randi(3); % pre target if no cues

% Choose target time for irregular
 if currentcondition==3 && ~qp && ~screeninfo.scopetest
    irrtartimeidx=alltrials(count,3);
    currirrtartimefr=timeinfo.irrtartimefr(irrtartimeidx);
elseif currentcondition==0 || qp || screeninfo.scopetest % Choose randomly for practice trials or scope tests
    irrtartimeidx=randi(5);
    currirrtartimefr=timeinfo.irrtartimefr(irrtartimeidx);
 end

% Choose gabor contrast
if ~screeninfo.scopetest
    %     if currentcondition==0 %visible practice
    %         currcontrast=1;
    %     else
    contrast=alltrials(count,1);
    currcontrast=visStiminfo.gcontrasts(contrast+1); % choose gabor contrast for this trial (+1 because 0 cannot index, so everything is shifted)
    %     end
end

if firstprac==1
    maskintensity=visStiminfo.pracmaskintensity;
    gaborpercent=visStiminfo.pracgaborpercent;
    noisetex=visStiminfo.pracnoisetex(randperm(length(visStiminfo.pracnoisetex)));
else
    maskintensity=visStiminfo.maskintensity;
    gaborpercent=visStiminfo.gaborpercent;
    noisetex=visStiminfo.noisetex(randperm(length(visStiminfo.noisetex)));
end

% Create Gabor embedded in mask
for i=1:timeinfo.TarDurFr
    if ~screeninfo.scopetest
        GaborAdjusted=currcontrast*visStiminfo.gaborcurr;
        randMask=(rand(visStiminfo.rectSize)-0.5)*2;
        targetStim=0.5+0.5*(maskintensity*randMask+GaborAdjusted*gaborpercent); % input: between -1 and 1, output between 0 and 1)
        visStiminfo.currTarget(i)=Screen('MakeTexture', screeninfo.win, targetStim);
    else
        noiseimg=0.01+maskintensity*(rand(visStiminfo.rectSize)-0.5);
        % make cue area of noise tex white
        xfrom=(length(visStiminfo.noiseimg)-visStiminfo.baseRect(3))/2;
        xto=xfrom+visStiminfo.baseRect(3);
        noiseimg(xfrom:xto,xfrom:xto)=visStiminfo.rectColorCue;
        visStiminfo.currTarget(i)=Screen('MakeTexture', screeninfo.win, noiseimg);
    end
end

% Choose target trigger
currtargtrig=[];

if currentcondition==1 % Rhythm
    irrISI=[0 0 0 0];
    currtargtrig=0;
    if trialorient==0 % Right
        currtargtrig=[currtargtrig 1];
    else
        currtargtrig=[currtargtrig 2];
    end
    if alltrials(count,1)==0 % Override for catch trials (don't follow trigger scheme)
        currtargtrig=[0 0 1];
    end
elseif currentcondition==3 % Irregular
    currtargtrig=1;
    if irrtartimeidx==1
        if trialorient==0 %right
            currtargtrig=[currtargtrig 0];
        else
            currtargtrig=[currtargtrig 5];
        end
    elseif irrtartimeidx==2
        if trialorient==0
            currtargtrig=[currtargtrig 1];
        else
            currtargtrig=[currtargtrig 6];
        end
    elseif irrtartimeidx==3
        if trialorient==0
            currtargtrig=[currtargtrig 2];
        else
            currtargtrig=[currtargtrig 7];
        end
    elseif irrtartimeidx==4
        if trialorient==0
            currtargtrig=[currtargtrig 3];
        else
            currtargtrig=[currtargtrig 8];
        end
    elseif irrtartimeidx==5
        if trialorient==0
            currtargtrig=[currtargtrig 4];
        else
            currtargtrig=[currtargtrig 9];
        end
    end
    if alltrials(count,1)==0 % Override for catch trials (don't follow trigger scheme)
        currtargtrig=[0 0 3];
    end
elseif currentcondition==2 % Interval
    irrISI=[0 0 0 0];
    currtargtrig=2;
    if trialorient==0
        currtargtrig=[currtargtrig 1];
    else
        currtargtrig=[currtargtrig 2];
    end
    if alltrials(count,1)==0 % Override for catch trials (don't follow trigger scheme)
        currtargtrig=[0 0 2];
    end

end

% Add contrast
if ~alltrials(count,1)==0 % If not a catch trial
   currtargtrig=[currtargtrig contrast];
end
currtargtrig=double(string(currtargtrig).join('')); % combine into one number

% Pick the rest of the triggers (use look up table):
% [onset cue1 cue2 cue3 warningsignal target question1 answereval11 answereval12 question2 answer21 answer22 answer23 answer24]
if currentcondition==1 && ~screeninfo.scopetest % Rhythm
    triggerlist=[31 41 51 61 71 currtargtrig 81 84 85 91 231 232 233 234];
elseif currentcondition==2 && ~screeninfo.scopetest % Interval
    triggerlist=[32 42 52 62 72 currtargtrig 82 86 87 92 241 242 243 244];
elseif currentcondition==3 && ~screeninfo.scopetest % Irregular
    triggerlist=[33 43 53 63 73 currtargtrig 83 88 89 93 250 251 252 253];
elseif screeninfo.scopetest
    triggerlist=[249 249 249 249 249 249 249 249 249 249 249 249 249 249];
else % for any condition without triggers
    triggerlist=zeros(14);
end


if currentcondition==2 % for interval choose a jittered between-interval time
    bintTime = randi([1400,1800])/1000; %%% use parameter %% % currently around 1600 to match total rhythm length
    timeinfo.bintTimefr=round(bintTime/screeninfo.ifi);
elseif currentcondition==3 || currentcondition==0 % determine all random ISI for irregular
    irrcondmet=0;
    while ~irrcondmet
        irrISI=timeinfo.irrwin(randperm(length(timeinfo.irrwin)));
        irrISI=irrISI(1:4);
        irrISIsort = sort(irrISI);
        if irrISIsort(2)-irrISIsort(1)>=timeinfo.mindiff && irrISIsort(3)-irrISIsort(2)>=timeinfo.mindiff && irrISIsort(4)-irrISIsort(1)>=timeinfo.mindiff
            irrcondmet=1;
        end
    end
    timeinfo.irrISIfr=round(irrISI/screeninfo.ifi);
end

% Show Noise

if ~justtar
    if currentcondition==1
        % Present Rhythm
        for fixtime=1:timeinfo.ITIFr(randi(length(timeinfo.ITIFr)))
            Screen('DrawTexture', screeninfo.win, noisetex(cnoise), [], [], 0);
            DrawFormattedText(screeninfo.win,blocklabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
            Screen('DrawLines', screeninfo.win, visStiminfo.fixCoords,...
                visStiminfo.lineWidthPix, screeninfo.black, [screeninfo.xCenter screeninfo.yCenter], 2);
            Screen('Flip', screeninfo.win);
            cnoise=cnoise+1;
            if fixtime==1
                io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, triggerlist(1));
            elseif fixtime==2
                io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, 0);
            end
        end
        for rhythms=1:3
            for CueDura=1:timeinfo.CueDurFr
                DrawFormattedText(screeninfo.win,blocklabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
                Screen('DrawTexture', screeninfo.win, visStiminfo.cuetex(ccue), [], [], 0);
                Screen('Flip', screeninfo.win);
                ccue=ccue+1;
                if CueDura==1
                    io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, triggerlist(1+rhythms));
                end
            end
            io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, 0);
            for rhythISI=1:timeinfo.rhyISIfr
                DrawFormattedText(screeninfo.win,blocklabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
                Screen('DrawTexture', screeninfo.win, noisetex(cnoise), [], [], 0);
                Screen('Flip', screeninfo.win);
                cnoise=cnoise+1;
            end
        end

        % Present Target Interval
        for WarnTime=1:timeinfo.CueDurFr
            DrawFormattedText(screeninfo.win,blocklabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
            Screen('DrawTexture', screeninfo.win, visStiminfo.warntex(cwarn), [], [], 0);
            Screen('Flip', screeninfo.win);
            cwarn=cwarn+1;
            if WarnTime==1
                io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, triggerlist(5));
            end
        end
        io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, 0);
        for rhythISI=1:timeinfo.rhyISIfr
            DrawFormattedText(screeninfo.win,blocklabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
            Screen('DrawTexture', screeninfo.win, noisetex(cnoise), [], [], 0);
            Screen('Flip', screeninfo.win);
            cnoise=cnoise+1;
        end

    elseif currentcondition==2
        for fixtime=1:timeinfo.ITIFr(randi(length(timeinfo.ITIFr)))
            Screen('DrawTexture', screeninfo.win, noisetex(cnoise), [], [], 0);
            DrawFormattedText(screeninfo.win,blocklabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
            Screen('DrawLines', screeninfo.win, visStiminfo.fixCoords,...
                visStiminfo.lineWidthPix, screeninfo.black, [screeninfo.xCenter screeninfo.yCenter], 2);
            Screen('Flip', screeninfo.win);
            cnoise=cnoise+1;
            if fixtime==1
                io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, triggerlist(1));
            elseif fixtime==2
                io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, 0);
            end
        end
        %     for intervals=1:2
        for CueDura=1:timeinfo.CueDurFr
            DrawFormattedText(screeninfo.win,blocklabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
            Screen('DrawTexture', screeninfo.win, visStiminfo.cuetex(ccue), [], [], 0);
            Screen('Flip', screeninfo.win);
            ccue=ccue+1;
            if CueDura==1
                 io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, triggerlist(2));
            end
        end
        io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, 0);
        for intervalISI=1:timeinfo.intISIfr
            DrawFormattedText(screeninfo.win,blocklabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
            Screen('DrawTexture', screeninfo.win, noisetex(cnoise), [], [], 0);
            Screen('Flip', screeninfo.win);
            cnoise=cnoise+1;
        end
        
        for CueDura=1:timeinfo.CueDurFr
            DrawFormattedText(screeninfo.win,blocklabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
            Screen('DrawTexture', screeninfo.win, visStiminfo.cuetex(ccue), [], [], 0);
            Screen('Flip', screeninfo.win);
            ccue=ccue+1;
            if CueDura==1
                io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, triggerlist(3));
            end
        end
        io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, 0);
        for intintISI=1:timeinfo.bintTimefr
            DrawFormattedText(screeninfo.win,blocklabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
            Screen('DrawTexture', screeninfo.win, noisetex(cnoise), [], [], 0);
            Screen('Flip', screeninfo.win);
            cnoise=cnoise+1;
        end

        %     end

        % Present Target Interval
        
        for WarnTime=1:timeinfo.CueDurFr
            DrawFormattedText(screeninfo.win,blocklabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
            Screen('DrawTexture', screeninfo.win, visStiminfo.warntex(cwarn), [], [], 0);
            Screen('Flip', screeninfo.win);
            cwarn=cwarn+1;
            if WarnTime==1
                io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, triggerlist(5));
            end
        end
        io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, 0);
        for intervalISI=1:timeinfo.intISIfr
            DrawFormattedText(screeninfo.win,blocklabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
            Screen('DrawTexture', screeninfo.win, noisetex(cnoise), [], [], 0);
            Screen('Flip', screeninfo.win);
            cnoise=cnoise+1;
        end

    else %for both condition 3 and practice condition 0
        for fixtime=1:timeinfo.ITIFr(randi(length(timeinfo.ITIFr)))
            Screen('DrawTexture', screeninfo.win,noisetex(cnoise), [], [], 0);
            DrawFormattedText(screeninfo.win,blocklabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
            Screen('DrawLines', screeninfo.win, visStiminfo.fixCoords,...
                visStiminfo.lineWidthPix, screeninfo.black, [screeninfo.xCenter screeninfo.yCenter], 2);
            Screen('Flip', screeninfo.win);
            cnoise=cnoise+1;
            if fixtime==1
                io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, triggerlist(1));
            elseif fixtime==2
                io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, 0);
            end
        end
    
        for irrint=1:3 % amount of intervals before target -1
           
            for CueDura=1:timeinfo.CueDurFr
                DrawFormattedText(screeninfo.win,blocklabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
                Screen('DrawTexture', screeninfo.win, visStiminfo.cuetex(ccue), [], [], 0);
                Screen('Flip', screeninfo.win);
                ccue=ccue+1;
                if CueDura==1
                    io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, triggerlist(1+irrint));
                end
            end
            io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, 0);
            for intervalISI=1:timeinfo.irrISIfr(irrint)
                DrawFormattedText(screeninfo.win,blocklabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
                Screen('DrawTexture', screeninfo.win, noisetex(cnoise), [], [], 0);
                Screen('Flip', screeninfo.win);
                cnoise=cnoise+1;
            end
        end

        % Present Target Interval
        for WarnTime=1:timeinfo.CueDurFr
            DrawFormattedText(screeninfo.win,blocklabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
            Screen('DrawTexture', screeninfo.win, visStiminfo.warntex(cwarn), [], [], 0);
            Screen('Flip', screeninfo.win);
            cwarn=cwarn+1;
            if WarnTime==1
                io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, triggerlist(5));
            end
        end
        io64(timeinfo.triggerPort, timeinfo.triggerPortAddress,0);
        for intervalISI=1:currirrtartimefr
            DrawFormattedText(screeninfo.win,blocklabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
            Screen('DrawTexture', screeninfo.win,noisetex(cnoise), [], [], 0);
            Screen('Flip', screeninfo.win);
            cnoise=cnoise+1;
        end
    end
else
    % Present Noise
    for pretartime=1:timeinfo.thrFr(pretarjitt)
        DrawFormattedText(screeninfo.win,blocklabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
        Screen('DrawTexture', screeninfo.win, noisetex(cnoise), [], [], 0);
        Screen('Flip', screeninfo.win);
        cnoise=cnoise+1;
    end
end

% Present Target
for TarTime=1:timeinfo.TarDurFr
    DrawFormattedText(screeninfo.win,blocklabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
    Screen('DrawTexture', screeninfo.win, visStiminfo.currTarget(TarTime));
    Screen('Flip', screeninfo.win);
    if TarTime==1
        io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, triggerlist(6));
    end
end

for QuestInt=1:timeinfo.TargQFr
    DrawFormattedText(screeninfo.win,blocklabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
    Screen('DrawTexture', screeninfo.win, noisetex(cnoise), [], [], 0);
    Screen('Flip', screeninfo.win);
    cnoise=cnoise+1;
end

DrawFormattedText(screeninfo.win,blocklabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
Screen('DrawLines', screeninfo.win, visStiminfo.fixCoords,...
    visStiminfo.lineWidthPix, screeninfo.black, [screeninfo.xCenter screeninfo.yCenter], 2);
Screen('Flip', screeninfo.win);

io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, 0);
pause(preqtime)

% Request Orientation Reponse
currtime=0;
warningobj=0;
subwarning=0;

if (alltrials(count,1)>0) && ~ screeninfo.speedrun && ~screeninfo.scopetest % skip questions for catch trials, speed runs and scope test
    startTime=GetSecs;
    orresponded=0;
    visresponded=0;
    io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, triggerlist(7));
    WaitSecs(0.01)
    io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, 0);
    while ~orresponded
        DrawFormattedText(screeninfo.win,blocklabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
        DrawFormattedText(screeninfo.win, textinfo.OrientRespPrompt, 'center', 'center', screeninfo.fontcolour);
        if (currtime-startTime) > timeinfo.maxresptime % if participant hasn't responded after max time, display warning
            DrawFormattedText(screeninfo.win,'Please answer.', 'center', screeninfo.yCenter+250, screeninfo.fontcolour);
            warningobj=1; % was a warning displayed during objective response?
        end
        Screen('Flip', screeninfo.win);
        [keyDown, respTime, keyName] = KbCheck;
        if keyDown
            if strcmp(KbName(keyName), textinfo.KeyLeft)==1
                RT=respTime-startTime;
                OrResp=1; % left
                orresponded=1;
                if trialorient==1 % left
                    io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, triggerlist(8));
                    OrRespEval=1;
                else
                    io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, triggerlist(9));
                    OrRespEval=0;
                end
            elseif strcmp(KbName(keyName), textinfo.KeyRight)==1
                RT=respTime-startTime;
                orresponded=1;
                OrResp=0; %Right
                if trialorient==0 %Right
                    io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, triggerlist(8));
                    OrRespEval=1;
                else
                    io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, triggerlist(9));
                    OrRespEval=0;
                end
            elseif strcmp(KbName(keyName),'e')==1
                DrawFormattedText(screeninfo.win,'Continue (c)? Exit (e)?', 'center', 'center', screeninfo.fontcolour);
                Screen('Flip', screeninfo.win);
                [~, keyNamethr, ~]=KbStrokeWait;
                if strcmp(KbName(keyNamethr),'e')
                    sca
                    return
                elseif strcmp(KbName(keyNamethr),'c')
                    OrResp=NaN; % if trial was paused in between, do not go back to record late answers but continue with next trial instead
                    OrRespEval=0;
                    VisResp=NaN;
                    visresponded=1;
                    RT=NaN;
                    DrawFormattedText(screeninfo.win,'Ready? Press any button.', 'center', 'center', screeninfo.fontcolour);
                    Screen('Flip', screeninfo.win);
                    KbStrokeWait
                    break
                end
            end
        end
        currtime=GetSecs;
    end
    WaitSecs(0.1)
    io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, 0);
    WaitSecs(0.1)
    % Give participant a speed reminder if they received a warning

    if warningobj
        DrawFormattedText(screeninfo.win,blocklabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
        DrawFormattedText(screeninfo.win,textinfo.warning, 'center', 'center', screeninfo.fontcolour);
        Screen('Flip', screeninfo.win);
        WaitSecs(1.5)
    end

    % Request Visibility Response
    if ~firstprac %for first visible practice, skip visiblity response
        DrawFormattedText(screeninfo.win,blocklabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
        Screen('Flip', screeninfo.win);
%         pause(0.2)
        currtime=0;
        startTime=GetSecs;
        io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, triggerlist(10));
        WaitSecs(0.01)
        io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, 0);
        while ~visresponded
            DrawFormattedText(screeninfo.win,blocklabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
            DrawFormattedText(screeninfo.win, textinfo.VisRespPrompt, 'center', 'center', screeninfo.fontcolour);
            if (currtime-startTime) > timeinfo.maxresptime % if participant hasn't responded after max time, display warning
                DrawFormattedText(screeninfo.win,'Please answer.', 'center', screeninfo.yCenter+250, screeninfo.fontcolour);
                subwarning=1; % was a warning displayed during objective response?
            end
            Screen('Flip', screeninfo.win);
            [keyDown, ~, keyName] = KbCheck;
            if keyDown
                if strcmp(KbName(keyName), textinfo.PAS1)==1
                    io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, triggerlist(11));
                    visresponded=1;
                    VisResp=0;
                elseif strcmp(KbName(keyName), textinfo.PAS2)==1
                    io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, triggerlist(12));
                    visresponded=1;
                    VisResp=1;
                elseif strcmp(KbName(keyName), textinfo.PAS3)==1
                    io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, triggerlist(13));
                    visresponded=1;
                    VisResp=2;
                elseif strcmp(KbName(keyName), textinfo.PAS4)==1
                    io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, triggerlist(14));
                    visresponded=1;
                    VisResp=3;
                elseif strcmp(KbName(keyName),'e')==1
                    DrawFormattedText(screeninfo.win,'Continue (c)? Exit (e)?', 'center', 'center', screeninfo.fontcolour);
                    Screen('Flip', screeninfo.win);
                    [~, keyNamethr, ~]=KbStrokeWait;
                    if strcmp(KbName(keyNamethr),'e')
                        sca
                        return
                    elseif strcmp(KbName(keyNamethr),'c')
                        VisResp=NaN; % don't record results for interrupted trials
                        RT=NaN;
                        DrawFormattedText(screeninfo.win,'Ready? Press any button.', 'center', 'center', screeninfo.fontcolour);
                        Screen('Flip', screeninfo.win);
                        KbStrokeWait
                        break
                    end
                end
            end
            currtime=GetSecs;
        end
    else
        VisResp=NaN;
    end
    WaitSecs(0.1)
    io64(timeinfo.triggerPort, timeinfo.triggerPortAddress, 0);
    % Give participant a speed reminder if they received a warning
    if subwarning
        DrawFormattedText(screeninfo.win,blocklabel,'center', screeninfo.yCenter-(0.25*screeninfo.axisy(2)),screeninfo.fontcolour);
        DrawFormattedText(screeninfo.win, textinfo.warning, 'center', 'center', screeninfo.fontcolour);
        Screen('Flip', screeninfo.win);
        WaitSecs(1.5)
    end
    Screen('Flip', screeninfo.win);
elseif screeninfo.speedrun % generate random answers for screeninfo.speedrun
    RT=0.4 + (2.4-0.4).*rand(1,1); % random reaction time between 0.4 and 2.4 seconds
    OrResp=randi(2)-1; % random orientation between 0 and 1
    if trialorient==1 % according evaluation
        OrRespEval=1;
    else
        OrRespEval=0;
    end
    VisResp=randi(4)-1;
else %for catch trials
    RT=NaN;
    OrResp=NaN;
    OrRespEval=NaN;
    VisResp=NaN;
    pause(timeinfo.ITI) %time in frames?
end
end

%% Adjust Intensity Function
% Insert current intensity, adjust and return new intensity.
% Increase/Decrease with arrows, confirm with enter
function [newintensity]=adjustintensity(screeninfo, gaborcontrast)
adjusting=1;
currlevel=gaborcontrast;
    while adjusting
    currentleveltext=sprintf('Current luminance level: %.2f',currlevel);
    DrawFormattedText(screeninfo.win,currentleveltext, 'center', 'center', screeninfo.fontcolour);
    Screen('Flip', screeninfo.win);
    
    [keyDown, ~, keyName] = KbCheck;
    if keyDown
        if strcmp(KbName(keyName), 'UpArrow')==1
            currlevel=currlevel+0.01;
            pause(0.05)
        elseif strcmp(KbName(keyName), 'DownArrow')==1
            currlevel=currlevel-0.01;
            pause(0.05)
        elseif strcmp(KbName(keyName), 'Return')==1 % confirm with enter
            newintensity=currlevel;
            adjusting=0;
        elseif strcmp(KbName(keyName),'e')==1 % exit with e
            sca
            return
        end
    
        if currlevel<0.01 % cannot go below 0.01
            currlevel=0.01;
        elseif currlevel>0.99 % or above 0.99
            currlevel=0.99;
        end
    end
    end
end