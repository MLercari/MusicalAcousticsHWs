%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Homework 4 exercise 4a
% Radiance pattern estimation directly from the signals
% We estimate the radiance pattern using the white noise source signal.
% The microphone signals are windowed accordingly to the first reflection
% time estimated in the previous exercises. 
%
% Musical Acoustic Course
% Mirco Pezzoli
% 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
close all
clc

%% Setup
fs = 48000;         % Sampling frequency
nfft = fs;          % Number of fft points
nMic = 24;          % Number of microphones
duration = 10;                               %[s] duration of sweep signal

R =   2.67;          % Distance between source and microphones

typeOfSignal = "noise"; % Noise
dir = strcat("Recordings/", typeOfSignal);% Signals directory

%% RADIANCE COMPUTATION THROUGH THE USE OF THE DIRECT RECORDED SIGNAL
% Window the signals according to the reflection time and estimate the
% radiance pattern. 

sig = zeros(fs*duration, nMic);   % Signal structure

% Early reflections attenuation. We consider only 2*ns+1 samples of the signal
ns = 150; % Ideal ns should be on the direct path
t = (0:1/fs:duration);
t = t(1:end-1);

TOA_directSignal = 0.0078;  % TOA
TOA_firstReflection =  0.011 ;% First reflection TOA

hann_win =  hann(2*ns +1);     % Window
w = zeros(length(t),1);
wind_size = (round(TOA_directSignal*fs) - ns):(round(TOA_directSignal*fs) + ns);
w(wind_size) = hann_win;

figure;
tiledlayout('flow');

x = audioread("input signals/noise.wav");   %input signal

for n = 1:nMic            % For each microphone signal
    % Load the signal
    y = audioread(strcat(dir, "/", num2str(n), ".wav"));
    y = y(1:duration*fs);
    
    % The window is applied to the signal from the 0 time instant.
    % Check the effect of differen window sizes.
    % with ns = 350 we include the reflection time instant, so we chose 150
    y_w = y.*w;
    
    % Plot the estimate impulse response
    [ir] = extractirnoise(x, y, nfft);
    plot(t, ir);
    
    nexttile
    % Plot the mic signal
    hold on
    plot(t,y) ;            % Plot the window over the signal
    plot(t, w*0.025);      % Scale the window for better visualize it
        
    % Plot the TOA using stem (see doc stem)
    y_toa1 = y(round(TOA_directSignal*fs));
    stem([TOA_directSignal], [y_toa1]);
    
    % Plot the first reflection TOA using stem (see doc stem)
    y_toa2 = y(round(TOA_firstReflection*fs));
    stem([TOA_firstReflection], [y_toa2]);
    
    % Plot the windowed signal with thicker line
    plot(t, y_w, 'LineWidth' , 0.5)
    hold off
    % Add a legend
    legend('signal', 'window', 'direct TOA', 'First reflection TOA', ...
        'windowed signal');
    xlim([0 0.05]);             % Limit the plot btw 0 and 0.05s
    xlabel('Time (sec)');
    title(['Mic: ', num2str(n)]);
    sig(:, n) = y_w; % Add current signal to the structure sig
end

%% Radiance estimation

SIG = fft(sig);% FFT of the windowed signal

rad_patt = % Compute the radiance pattern magnitude

%% Radiance pattern plot

% Frequencies must be centered at a frequency bins
ctr_freqs = [150 250 500 750 1000 1500 2000 5000 7500 10000];

% Microphone angle wrt the speaker frontal face
angs = 

% Plot the estimated radiance pattern using the provided function
% radianceplot

radianceplot(ctr_freqs, rad_patt, angs, 'noise direct: ');

