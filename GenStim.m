function [Audio,start_times] = GenStim(freq,phase,trial_num,tone_num,dev_type,deviations,dev_tone,fs)

% Set a timer to continually update with the times of the different tones
% Start at 0
    % Add durations to the timer
% Add the start times to a separate array
% If the start time belongs to a deviant, add it to another array


%%

durTone = 0.40; % duration of TONE stimulus (seconds)
durDev = 0.20; % duration of the timing deviation
durQuiet = 0.50;

tTone = (1/fs:1/fs:durTone)'; % time vector for each individual stim
tDev = (1/fs:1/fs:durDev)'; % time vector for a time deviation

tones = [];
pauses = [];
pattern = []; % combination of the three tones [tone, tone, tone]
Audio = []; % complete audio file [pat, pat, pat, ..., pat]
freq_tracker = []; %Keeps track of which frequencies are being used
dev_freqs = [];
demo_num = 2; % number of demonstration patterns to establish the tones before deviants are introduced

dev_times = [];

%% Generate the tones

for f = 1:(tone_num+1) % Based on the number of tones in the pattern, last tone is the deviant
    
    tone = sin(2*pi*freq(f)*tTone+phase); % Define the tone
    onoffRamp = 0.016; % duration of cosine-squared ramps (seconds)
    winTone = tukeywin(length(tone),2*onoffRamp/durTone);
    tone = tone.*winTone; % pure tone stimulus
    
    % Concatenate the tones together into a pattern
    tones(f,:) = tone;
    norm_freqs(1:tone_num) = freq(1:tone_num);
    % Define pauses between tones
    if f > 1 && f <= tone_num
        durPause = randi([0 3])/10; % random pause between tones
        pauses(f-1) = durPause;
    end
end
%% Combine tones into a pattern

durPattern = durTone*tone_num+durQuiet; % duration of the pattern (quiet included)
% tP = (0:1/fs:durPattern)'; % time vector for the pattern

pattern_norm = tones(1,:); % Begin the definition of the pattern by inserting the first tone
pattern_times = [0];

for ii = 2:tone_num % Add the inter-tone pauses and subsequent tones to the pattern
    tPause = (1/fs:1/fs:pauses(ii-1)); % time vector for the pause between tones
    Pause = zeros(length(tPause),1)';
    pattern_norm = [pattern_norm,Pause,tones(ii,:)];
    ToneTime = pattern_times(end)+durTone+pauses(ii-1);
    pattern_times = [pattern_times,ToneTime];
end

%% Define deviant patterns

if dev_type == 1 % for timing deviations:
    
    pattern_dev = tones(1,:); % the tones themselves will remain the same
    dev_pat_times = 0;
    
    for jj = 2:tone_num % modifies the inter-tone pauses
        
        if jj == 2
            tPause = (1/fs:1/fs:(pauses(jj-1)+durDev));
            ToneTime = dev_pat_times(end)+durTone+(pauses(jj-1)+durDev);
        else
            tPause = (1/fs:1/fs:pauses(jj-1));
            ToneTime = dev_pat_times(end)+durTone+pauses(jj-1);
        end
        
        Pause = zeros(length(tPause),1)';
        pattern_dev = [pattern_dev,Pause,tones(jj,:)];
        dev_pat_times = [dev_pat_times,ToneTime];
    end
        
elseif dev_type == 2 % for frequency deviations:
    
    if dev_tone == 1 %Define the first tone
        pattern_dev = tones(tone_num+1,:);
        dev_freqs(1) = freq(tone_num+1);
    else
        pattern_dev = tones(1,:);
        dev_freqs(1) = freq(1);
    end
        
    for jj = 2:tone_num
        tPause = (1/fs:1/fs:(pauses(jj-1)));
        Pause = zeros(length(tPause),1)';
        if jj == dev_tone
            pattern_dev = [pattern_dev,Pause,tones(tone_num+1,:)];
            dev_freqs(jj) = freq(tone_num+1);
        else
            pattern_dev = [pattern_dev,Pause,tones(jj,:)];
            dev_freqs(jj) = freq(jj);
        end        
    end
end

 
%% Combine patterns into one stimulus

deviations = deviations+demo_num; %This will bump the deviation positions back so that the first two demos will be the normal pattern

for g = 1:trial_num+demo_num
    
    durQuiet = 0.50; % duration of quiet between patterns
    
    if ismember(g,deviations) && dev_type == 1 % If there is a time deviation, adjust the timing
        durQuiet = durQuiet-durDev; %shorten the time between patterns
        tQ = (0:1/fs:durQuiet)'; % time vector for the quiet
        pattern = pattern_dev;
        freq_tracker = [freq_tracker;norm_freqs(1:tone_num)];
        dev_times = [dev_times];
        pat_timer = dev_pat_times;
    elseif ismember(g,deviations) && dev_type == 2 % If this is a frequency deviation
        tQ = (0:1/fs:durQuiet)'; % time vector for the quiet
        pattern = pattern_norm;
        freq_tracker = [freq_tracker;dev_freqs];
        dev_times = [];
        pat_timer = pattern_times;
    else % otherwise, combine the tones into a normal pattern
        tQ = (0:1/fs:durQuiet)'; % time vector for the quiet
        pattern = pattern_norm;
        freq_tracker = [freq_tracker;norm_freqs(1:tone_num)];
        pat_timer = pattern_times;
    end
    
    freq_tracker = [freq_tracker;norm_freqs(1:tone_num)];
    Quiet = zeros(length(tQ),1);
    
    if g == 1 %Add the patterns together
        Audio = pattern;
        start_times = pattern_times;
    else
        Audio = [Audio, Quiet', pattern];
        current_pattern_time = start_times(end)+durTone+durQuiet+pat_timer;
        start_times = [start_times,current_pattern_time];
    end
    
    if ismember(g,deviations) && dev_type == 1
        dev_times = [dev_times,current_pattern_time(1)];
    elseif ismember(g,deviations) && dev_type == 2 % If there is a deviation, add it to the time tracker
        dev_times = [dev_times,current_pattern_time(dev_tone)];
    end
 
end

start_times(1) = 1/fs;
%% Normalize signal

% Audio = zeros(size(durStim));
% Audio(k0:k0+length(tone)-1) = tone; % create TONE vector of same size as NOISE

% Normalize to ±1
% Audio = Audio/max(max(abs([Audio,NOISE])));

%% Plot signal

time_vec = [1/fs:1/fs:length(Audio)/fs];
t_locs = round(start_times*44100);
stim_timer = zeros(1,length(Audio));
stim_timer(t_locs) = 1;

d_locs = round(dev_times*44100);
dev_timer = zeros(1,length(Audio));
dev_timer(d_locs) = 1;

plot(time_vec,Audio,'r')
hold on
% plot(time_vec,stim_timer,'k')
plot(time_vec,dev_timer,'k','LineWidth',2)
grid on
grid minor


end
