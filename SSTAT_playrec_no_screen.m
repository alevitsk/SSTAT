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

siglvl=60;
N = 3; % number of audio stimuli
trial_num = 20;
devs = randperm(trial_num,10); % Define the positions of two deviations
dev_pos = randi(N);
audio_level = 60;

    %% Generate stimulus audio
    fprintf('Generating stimulus audio...');
    phase=0;
    stim_order = randi(7,N,1);
    
    freq_dev = randi([1,2]); % Decides whether to increase or lower the frequency of the deviant tone
    if stim_order(dev_pos) < 2
        stim_order(N+1) = stim_order(dev_pos)+1;
    elseif stim_order(dev_pos) == length(freq_list) 
        stim_order(N+1) = stim_order(dev_pos)-1;        
    elseif freq_dev == 1
        stim_order(N+1) = stim_order(dev_pos)+1;
    elseif freq_dev == 2
        stim_order(N+1) = stim_order(dev_pos)-1;        
    end
    
    stim_freq = freq_list(stim_order); %Defining which frequencies are used for the tones
    [xt] = GenStim(stim_freq,phase,trial_num,N,dev_type,devs,dev_pos,fs);
    amp = db2mag(audio_level-85);  %Need to confirm attenuation levels for the speakers
    fprintf('DONE\n');
    
    %% Experimental trial mock up
%     [RT] = PresentStim(xt,nt,amp,t0,win,scrY,deviceIndex,tt);
    
    STIM = [amp*xt];
    playrec('play',STIM',1);

    wavwrite('SSTAT_Audio.wav',STIM,fs);
