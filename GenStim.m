function [Audio] = GenStim(freq,phase,trial_num,tone_num,dev_type,deviations,dev_tone,fs)

% First, initialize the tones and the quiet
% Set a loop to create the tones
% Have the stims created twice, and play back-to-back
% Set a loop to create the pattern based on the number of individual tones that you will be hearing ('N' value).
% Create a deviation, or don't
% Set the timing of the deviation to .05 seconds
% In a loop, put all of them together

%What's the best way to do this?

% Combine the stimuli into one pattern.
% Combine the pattern

%%

% if nargin<8
%     fs = 44100;
% end

durStim = 0.40; % duration of TONE stimulus (seconds)
durTotal = 45.0; %duration of the entire stimulus
durDev = 0.050; % duration of the timing deviation

tTone = (0:1/fs:durStim)'; % time vector for each individual stim
tTotal = (0:1/fs:durTotal)'; % time vector for the complete stream
tDev = (0:1/fs:durDev)'; % time vector for a time deviation

tones = [];
pauses = [];
pattern = []; % combination of the three tones [tone, tone, tone]
Audio = []; % complete audio file [pat, pat, pat, ..., pat]
freq_tracker = []; %Keeps track of which frequencies are being used
dev_freqs = [];

%% Generate the tones

for f = 1:(tone_num+1) % Based on the number of tones in the pattern, last tone is the deviant
    
    tone = sin(2*pi*freq(f)*tTone+phase); % Define the tone
    
    onoffRamp = 0.016; % duration of cosine-squared ramps (seconds)
    winTone = tukeywin(length(tone),2*onoffRamp/durStim);
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

% durPattern = durStim*tone_num+durQuiet; % duration of the pattern (quiet included)
% tP = (0:1/fs:durPattern)'; % time vector for the pattern
% 
% start_times = (0:durPattern:durTotal);

pattern_norm = tones(1,:); % Begin the definition of the pattern by inserting the first tone

for ii = 2:tone_num % Add the inter-tone pauses and subsequent tones to the pattern
    tPause = (0:1/fs:pauses(ii-1)); % time vector for the pause between tones
    Pause = zeros(length(tPause),1)';
    pattern_norm = [pattern_norm,Pause,tones(ii,:)];
end

%% Define deviant patterns

if dev_type == 1 % for timing deviations:
    
    pattern_dev = tones(1,:); % the tones themselves will remain the same
    
    for jj = 2:tone_num % modifies the inter-tone pauses
        
        if jj ==2
            tPause = (0:1/fs:(pauses(jj-1)+tDev));    
        end
        
        tPause = (0:1/fs:pauses(ii-1)); 
        Pause = zeros(length(tPause),1)';
        pattern_dev = [pattern_dev,Pause,tones(jj,:)];
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
        tPause = (0:1/fs:(pauses(jj-1)));
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

for g = 1:trial_num
    
    durQuiet = 0.50; % duration of quiet between patterns
    
    if ismember(g,deviations) && dev_type == 1 % If there is a time deviation, adjust the timing
        durQuiet = durQuiet-durDev; %shorten the time between patterns
        tQ = (0:1/fs:durQuiet)'; % time vector for the quiet
        pattern = pattern_dev;
        freq_tracker = [freq_tracker;norm_freqs(1:tone_num)];
    elseif ismember(g,deviations) && dev_type == 2
        tQ = (0:1/fs:durQuiet)'; % time vector for the quiet
        pattern = pattern_norm;
        freq_tracker = [freq_tracker;dev_freqs];
    else % otherwise, combine the tones into a normal pattern
        tQ = (0:1/fs:durQuiet)'; % time vector for the quiet
        pattern = pattern_dev;
        freq_tracker = [freq_tracker;norm_freqs(1:tone_num)];
    end
    
    freq_tracker = [freq_tracker;norm_freqs(1:tone_num)];
    
    Quiet = zeros(length(tQ),1);
    
    if g == 1 %Add the patterns together
        Audio = pattern;
    else
        Audio = [Audio, Quiet', pattern];
    end
end


%% Normalize signal

% Audio = zeros(size(durStim));
% Audio(k0:k0+length(tone)-1) = tone; % create TONE vector of same size as NOISE

% Normalize to ±1
% Audio = Audio/max(max(abs([Audio,NOISE])));
end
