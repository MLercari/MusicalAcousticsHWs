%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Homework 4 exercise 3b
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

nMic = 24;
typeOfSignal = "sweep";       % Sweep
dir = "Recordings/sweep/";  % File directory

fs = 48000;                                 % Sampling frequency
c = 343.8;                     %[m]/[s]            
duration = 10;                               %[s] duration of sweep signal

% This is the input signal. The following function computes 1/FFT(sweep). 
% We will thus need to MULTIPLY this signal for FFT(output) in order to 
% apply deconvolution (i.e. computing the radiation pattern)
% Use the provided synthSweep function.
[sweep, invsweepfft ,sweepRate] = synthSweep(duration,fs,50,22e3,0);
figure;
tiledlayout('flow');

direct = zeros(1, nMic); %store the direct path length here
reflect = zeros(1, nMic); %store the first reflection length here

for n = 1:nMic    % For each microphone signal
    
    % Load the signal
     x = audioread(strcat(dir, num2str(n), '.wav'));
     x = x(1:duration*fs);
     
    % Compute the impulse response using the function extractirsweep
    [ir] = extractirsweep(x, invsweepfft);

    % Setting up time scale for computed ir
    t =  (0:1/fs:length(ir)/fs);
    t = t(1:end -1);
    
    % Find the first impulse of the impulse response around
    % the peak position (0.0077 s)
    [pks, locs] = findpeaks(ir(0.005*fs:0.01*fs), t(0.005*fs:0.01*fs) , 'NPeaks', 1 , 'SortStr' , ...
     'descend' , 'MinPeakHeight' , 0);
    
    % Double check if we are finding back the correct direct path length
    directPathTimeOfArrival = locs;
    directPathLength = directPathTimeOfArrival*c;
    
    % Find the second peak of the impulse response around
    % the peak position (0.01 s)
    [pks2, locs2] = findpeaks(ir(0.009*fs:0.012*fs), t(0.009*fs:0.012*fs) , 'NPeaks', 1 , 'SortStr' , ...
        'descend' , 'MinPeakHeight' , 0);
    
    firstReflectionTimeOfArrival = locs2; %[s]
    firstReflectionPathLength = firstReflectionTimeOfArrival*c - directPathLength; %[m]
    
    % Plot the estimated impulse response
    nexttile
    hold on
    plot(t, ir);
    y_stem = NaN(1,length(t));
    y_stem( round(locs*fs) ) = ir( round(locs*fs) );
    stem( t , y_stem)
    
    y_stem2 = NaN(1,length(t));
    y_stem2( round(locs2*fs) ) = ir( round(locs2*fs) );
    stem( t , y_stem2)
    
    hold off
    xlim([0 0.05]);
    xlabel('Time (sec)');
    title(['Mic: ', num2str(n)]);
    
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
% 
 fprintf(sprintf('Direct path length %f m\n', meanDirect));

%% MIC TO WALLS DISTANCE COMPUTATION (ROOM SIZE)

% Put here the difference between first reflection and direct sound time of
% arrivals (TOA). By visual inspection we find that the second peak in
% the autocorrelation is around 0.010983 s
  
% Compute the distance from the reflectors
  distance = mean(reflect); 

% Print on screen the estimated distance from the walls
 fprintf(sprintf('Average distance between first path and reflector %f m\n', ...
     distance));
