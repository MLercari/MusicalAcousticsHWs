%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Homework 4 exercise 3a
% Evaluation of reflections in the signals
% We compute the impulse response using the white noise source signal
% from which the first direct path is computed. In addition by a visual 
% inspection we estimate the first reflection time of arrival. Both the 
% distance from the source and the distance from the reflectors are 
% computed using the estimated times.
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
addpath('input signals')

nMic = 1;              % Number of microphones
c = 343.8; % [m]/[s]

typeOfSignal = "noise"; % Noise
dir = "Recordings/noise/" ; % File directory

inputSignalDir = "input signals/" ;            % Source signal directory
inputSignalFileName = "noise.wav" ;    % Source signal name

%% Impulse response estimation

% The input signal must be known, load the input signal

[x, fs] = audioread(strcat(inputSignalDir, inputSignalFileName));
nfft = fs;              % Number of fft points
t = (0:1/fs:nfft/fs);   % Time axis
t = t(1:end-1);    

figure;
tiledlayout('flow');
for n = 1:nMic            % For each microphone signal
    
    % Load the signal
    y = audioread( strcat( dir, num2str(n), ".wav"));
    
    % Cut the recordings accoding to input signal length
    y = y(1:length(x));
    
    % Compute the impulse response using the function extractirnoise
    [ir] = extractirnoise(x, y, nfft);
    
    % Find the first (and highest) impulse of the impulse response
    [pks, locs] = findpeaks(ir, t , 'NPeaks', 1 , 'SortStr' , ...
        'descend' , 'MinPeakHeight' , 0);
    
    % Double check if we are finding back the correct direct path length
    
    directPathTimeOfArrival = locs; %[s]
    directPathLength = directPathTimeOfArrival*c;
    
    % Plot the estimate impulse response
    nexttile
    hold on
    plot(t, ir);
    y_stem = NaN(1,length(t));
    y_stem( round(locs*fs) ) = ir( round(locs*fs) );
    stem( t , y_stem)
    xlim([0 0.05]);
    xlabel('Time (sec)');
    title(['Mic: ',num2str(n)]);
    
end

%% SPEAKER TO MIC DISTANCE (DIRECT PATH LENGTH)
% Print on screen the estimated distance from the source

fprintf(sprintf('Direct path length %f m\n', directPathLength));

%% MIC TO REFLECTORS DISTANCE COMPUTATION

% Put here the difference between first reflection and direct sound time of
% arrivals (TOA)

% Inspecting the impulse responses determine the delay of the first
% reflection
  delay = 0.010983 - directPathTimeOfArrival; %[s]
  
% Compute the distance from the reflectors
  distance = delay*c; 

% Print on screen the estimated distance from the reflectors
fprintf(sprintf('Average distance between first path and reflector %f m\n', distance));
