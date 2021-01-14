%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Homework 4 exercise 4b
% Radiance pattern estimation directly from the signals
% We estimate the radiance pattern using the sine sweep source signal.
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
nMic = 1;          % Number of microphones
duration = 10;                               %[s] duration of sweep signal
c = 343.8;                     %[m]/[s]      
R =  2.67;          % Distance between source and microphones

typeOfSignal = "sweep";      % Sweep

dir = 'Recordings/sweep/';      % Signals directory

%% RADIANCE COMPUTATION THROUGH THE USE OF THE DIRECT RECORDED SIGNAL
% Window the signals according to the reflection time and estimate the
% radiance pattern. 
sig = [];   % Signal structure

% Early reflections attenuation. We consider only 2*ns+1 samples of the signal
%ns = 350;
ns = 300; % Ideal ns should be on the direct path
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

for n = 1:nMic            % For each microphone signal
    % Load the signal
    
    y = audioread(strcat(dir, num2str(n), '.wav'));
    y = y(1:duration*fs);
    % The window is applied to the signal from the 0 time instant.
    % Check the effect of differen window sizes.
    w = w.*0.025;
    y_w = y.*w;
    
    nexttile    % Plot the estimate impulse response
    [sweep, invsweepfft ,sweepRate] = synthSweep(duration,fs,50,22e3,0);
    [ir] = extractirsweep(y_w, invsweepfft);
    time = t(1:length(ir));
    plot(time, ir);
    xlim([0 0.5])
    
    nexttile
    % Plot the mic signal
    hold on
    plot(t,y)             % Plot the window over the signal
    plot(t, w)
    % Plot the TOA using stem (see doc stem)
    x_stem = y;
    y_stem1 = zeros(length(t),1);
    y_stem1(round(TOA_directSignal*fs)) = ir(round(TOA_directSignal*fs));
    % Plot the first reflection TOA using stem (see doc stem)
    y_stem2 = zeros(length(t),1);
    y_stem2(round(TOA_firstReflection*fs)) = ir(round(TOA_firstReflection*fs));
    % Plot the windowed signal with thicker line
    plot(t, y_w, 'LineWidth' , 0.5)
    hold off
    % Add a legend
    
    xlim([0 0.05]);             % Limit the plot btw 0 and 0.05s
    xlabel('Time (sec)');
    
    sig =  y_w;           % Add current signal to the structure sig
end
%% Radiance estimation
%{
SIG = % FFT of the windowed signal

rad_patt = % Compute the radiance pattern magnitude

%% Radiance pattern plot
% Frequencies must be centered at a frequency bins
ctr_freqs = [150 250 500 750 1000 1500 2000 5000 7500 10000];

% Microphone angle wrt the speaker frontal face
angs = 

% Plot the estimated radiance pattern using the provided function
% radianceplot
radianceplot(ctr_freqs, rad_patt, angs, 'sweep direct: '); 
%}
