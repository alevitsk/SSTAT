clear all; close all; clc;
commandwindow;

%% Initialize sound card
% Step One: Connect to and properly initialize sound card
ptb = genpath('/home/andrew/MATLAB/Psychtoolbox/');
pr = genpath('/home/andrew/MATLAB/PlayrecForMatlab-playrec-7c15bf6/');
addpath(ptb);
addpath(pr);
fprintf('Initializing connection to sound card...\n')

% %% Initialize sound card (RME)
% % Step One: Connect to and properly initialize RME sound card
fprintf('Initializing connection to sound card...\n')
Devices=playrec('getDevices');
if isempty(Devices)
    error(sprintf('There are no devices available using the selected host APIs.\nPlease make sure the RME is powered on!')); %#ok<SPERR>
else
    i=1;
    while ~strcmp(Devices(i).name,'default') && i <= length(Devices)
        i=i+1;
    end
end
stimchanList=[1];
fs = Devices(i).defaultSampleRate;
playDev = Devices(i).deviceID;
playrec('init',fs,playDev,-1,32,-1);
fprintf('Success! Connected to %s.\n', Devices(i).name);

todayStr = datestr(now,'yyyymmdd');
inputInfo=inputdlg({'subject ID: '});
sID = inputInfo{1};

deviationflag = 1; %Define whether the deviant tones will be time-displaced or have different frequencies
while deviationflag == 1
    dev_type = input('Please enter deviation type (F or T for frequency/time, respectively:', 's');
    switch dev_type
        case {'T', 't', 'time', 'Time', 'TIME'}
            dev_type = 1;
            deviationflag = 0;
        case {'F', 'f', 'freq', 'frequency', 'Frequency', 'FREQUENCY'}
            dev_type = 2;
            deviationflag = 0;
        otherwise
            fprintf(2, 'Unrecognized answer! Try again!');
    end
end

%% Define experimental parameters
freq_list = [500,750,1000,1250,1500,1750,2000];

N = 3; % number of audio stimuli
trial_num = 20; % Number of pattern repetitions
devs = randperm(trial_num,10); % Define the positions of two deviations
dev_pos = randi(N); % Defines which tone will be deviant
audio_level = 60;

%% Generate stimulus audio
fprintf('Generating stimulus audio...');
phase=0;
stim_order = randi(7,N,1); % Decides which frequencies will be used for the tones (Nx1 frequencies from the list of 7)

% Adds an extra item to define the deviant tone if we have a
% frequency deviation
% 'If' statements decide whether to increase or lower the frequency of the
% deviant tone

freq_dev = randi([1,2]);
if stim_order(dev_pos) < 2 % If it's too early in the list
    stim_order(N+1) = stim_order(dev_pos)+1; % Moves the frequency farther along the list (higher frequency)
elseif stim_order(dev_pos) == length(freq_list) % If it's too far in the list
    stim_order(N+1) = stim_order(dev_pos)-1; % Moves the frequency lower in the list (lower frequency)
elseif freq_dev == 1
    stim_order(N+1) = stim_order(dev_pos)+1;
elseif freq_dev == 2
    stim_order(N+1) = stim_order(dev_pos)-1;
end

stim_freq = freq_list(stim_order); %Defining which frequencies are used for the tones
[xt,t] = GenStim(stim_freq,phase,trial_num,N,dev_type,devs,dev_pos,fs);
amp = db2mag(audio_level-85);  %Need to confirm attenuation levels for the speakers
fprintf('DONE\n');

%% Experimental trial mock up
%     [RT] = PresentStim(xt,nt,amp,t0,win,scrY,deviceIndex,tt);

STIM = [amp*xt];
playrec('play',STIM',1);

audiowrite('SSTAT_Audio.wav',STIM,fs);

%% Plot the audio

time_vec = [1/fs:1/fs:length(STIM)/fs];
t_locs = round(t*44100);
stim_timer = zeros(1,length(STIM));
stim_timer(t_locs) = 1;

plot(time_vec,xt,'g')
hold on
plot(time_vec,stim_timer,'k')
grid on
grid minor

%%  Scoring

Screen('FillRect',win,black);
Screen('Flip',win);

KbQueueRelease(deviceIndex);
WaitSecs(2);
ShowCursor;
ListenChar(0);
sca;

save(['Results/', todayStr '_' sID '_SRFMData.mat'],'SRFMData');

t=toc;
fprintf('Experiment complete. Total time elapsed %d min %d sec.\n',floor(t/60),floor(mod(t,60)));
