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

addpath('Functions')
addpath('Recordings')
addpath('input signals')


fs = 48000;         % Sampling frequency
nfft = fs;          % Number of fft points
nMic = 24;          % Number of microphones
duration = 3;                               %[s] duration of sweep signal
c = 343.8;                     %[m]/[s]
R =   2.67;          % Distance between source and microphones
typeOfSignal = "noise"; % Noise
dir = strcat("Recordings/", typeOfSignal , "/" ); % Signals directory

%% RADIANCE COMPUTATION THROUGH THE USE OF THE DIRECT RECORDED SIGNAL
% Window the signals according to the reflection time and estimate the
% radiance pattern. 

sig = zeros(fs*duration, nMic);   % Signal structure

% Early reflections attenuation. We consider only 2*ns+1 samples of the signal
ns = 350; % with this ns the first reflection time is included
%ns = 120; %with this ns the first reflection time is exluded BUT the
%frequency resolution is too low!!

%time vector
t = 0:1/fs:duration ;
t = t(1:end-1);

%from previous excercise 3a we know that the mean values are

TOA_directSignal = 0.008;  % direct signal TOA
TOA_firstReflection =  0.008 + 0.0032; % First reflection TOA

hann_win =  hann(2*ns +1);     % Window
w = zeros(length(t),1);
wind_size = (round(TOA_directSignal*fs) - ns):(round(TOA_directSignal*fs) + ns);
w(wind_size) = hann_win;

figure;
tiledlayout('flow');

%source signal
x = audioread("input signals/noise.wav");   %input signal

for n = 1:nMic            % For each microphone signal
    
    % Load the signal (which is longer than x)
    y = audioread(strcat(dir, "/", num2str(n), ".wav"));
    y = y(1:duration*fs);
    
    % The window is applied to the signal from the 0 time instant (only if ns = 350).
    % Check the effect of differen window sizes.
    % with ns = 350 we include the reflection time instant, so we chose 150
    
    y_w = y.*w;
    
    nexttile
    % Plot the mic signal
    hold on
    plot(t,y) ;            % Plot the window over the signal
    plot(t, w*0.025);      % Scale the window for better visualize it
        
    % Plot the TOA using stem (see doc stem)
    x_stem = t;
    y_stem1 = NaN(length(t),1);
    y_stem1(round(TOA_directSignal*fs)) = y(round(TOA_directSignal*fs));
    stem(x_stem, y_stem1 );
    
    % Plot the first reflection TOA using stem (see doc stem)
    y_stem2 = NaN(length(t),1);
    y_stem2(round(TOA_firstReflection*fs)) = y(round(TOA_firstReflection*fs));
    stem(x_stem, y_stem2 );
    
    % Plot the windowed signal with thicker line
    plot(t, y_w, 'LineWidth' , 0.5)
    hold off
    
    xlim([0 0.02]);             % Limit the plot btw 0 and 0.05s
    xlabel('Time (sec)');
    title(['Mic: ', num2str(n)]);
    sig(:, n) = y_w; % Add current signal to the structure sig
    
end
    % Add a legend
    lgd = legend('signal', 'window', 'direct TOA', 'First reflection TOA', ...
        'windowed signal');

    lgd.Layout.Tile = "south";
    
%% Radiance estimation

SIG = fft(sig, nfft); % FFT of the windowed signal (important to set nfft)

S = fft(x , nfft); % input fft: important to set nfft


f = 0:fs; 
f = f(1:end -1);
omega = f.*(2*pi); %coordinate for the green function - use this? or f?

G = exp(-1i*(f./c)*R)./(4*pi*R); %green function 
rad_patt = zeros(fs, nMic);

for i = 1:nMic
    
rad_patt(:,i) = SIG(:,i)./(S.*G'); % Compute the radiance pattern magnitude

end


%% Radiance pattern plot
% Frequencies must be centered at a frequency bins
ctr_freqs = [150 250 500 750 1000 1500 2000 5000 7500 10000];

% Microphone angle wrt the speaker frontal face
angs = zeros(1,nMic);

for ii = 2:nMic
    angs(1,ii) = angs(1,ii-1) + 15;
end

% Plot the estimated radiance pattern using the provided function
% radianceplot
radianceplot(ctr_freqs, rad_patt, angs, 'noise direct: '); 

