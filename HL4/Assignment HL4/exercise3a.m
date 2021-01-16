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

nMic = 24;              % Number of microphones
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

direct = zeros(1, nMic); %store the direct path length here
reflect = zeros(1, nMic); %store the first reflection length here
for n = 1:nMic            % For each microphone signal
    
    % Load the signal
    y = audioread( strcat( dir, num2str(n), ".wav"));
    
    % Cut the recordings accoding to input signal length
    y = y(1:length(x));
    
    % Compute the impulse response using the function extractirnoise
    [ir] = extractirnoise(x, y, nfft);
    
    % Find the first (and highest) impulse of the impulse response around
    % the peak position (0.0077 s)
    [pks, locs] = findpeaks(ir(0.005*fs:0.01*fs), t(0.005*fs:0.01*fs) , 'NPeaks', 1 , 'SortStr' , ...
        'descend' , 'MinPeakHeight' , 0);
    
    % Double check if we are finding back the correct direct path length
    
    directPathTimeOfArrival = locs; %[s]
    
    directPathLength = directPathTimeOfArrival*c;
    
    [pks2, locs2] = findpeaks(ir(0.009*fs:0.012*fs), t(0.009*fs:0.012*fs) , 'NPeaks', 1 , 'SortStr' , ...
        'descend' , 'MinPeakHeight' , 0);
    
    firstReflectionTimeOfArrival = locs2; %[s]
    firstReflectionPathLength = firstReflectionTimeOfArrival*c - directPathLength; %[m]
    
    % Plot the estimate impulse response
    nexttile
    hold on
    plot(t, ir);
    y_stem = NaN(1,length(t));
    y_stem( round(locs*fs) ) = ir( round(locs*fs) );
    stem( t , y_stem)
    
    y_stem2 = NaN(1,length(t));
    y_stem2( round(locs2*fs) ) = ir( round(locs2*fs) );
    stem( t , y_stem2)
    
    xlim([0 0.05]);
    xlabel('Time (sec)');
    title(['Mic: ',num2str(n)]);
    
    direct(n) = directPathLength;
    reflect(n) = firstReflectionPathLength;
end

%% SPEAKER TO MIC DISTANCE (DIRECT PATH LENGTH)
% Print on screen the estimated distance from the source

meanDirect = mean(direct); %the mean value of path length

 figure(2);
 stem(1:nMic, direct )
 ylim([2.5 3])
 xlabel('Measurement'), ylabel('Distance highest peak')
 
fprintf(sprintf('Direct path length %f m\n', meanDirect));

%% MIC TO REFLECTORS DISTANCE COMPUTATION

% Put here the difference between first reflection and direct sound time of
% arrivals (TOA). By visual inspection we find that the second peak in
% the autocorrelation is around 0.010983 s
  
% Compute the distance from the reflectors
  distance = mean(reflect); 

% Print on screen the estimated distance from the reflectors
 fprintf(sprintf('Average distance between first path and reflector %f m\n', distance));
