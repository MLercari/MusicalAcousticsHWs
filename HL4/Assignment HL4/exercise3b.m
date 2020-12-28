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

nMic = 
typeOfSignal =                     % Sweep
dir =   % File directory

fs = 48000;                                 % Sampling frequency
speed_of_sound = 343.8;                     %[m]/[s]            
duration = 10;                               %[s] duration of sweep signal

% This is the input signal. The following function computes 1/FFT(sweep). 
% We will thus need to MULTIPLY this signal for FFT(output) in order to 
% apply deconvolution (i.e. computing the radiation pattern)
% Use the provided synthSweep function.
figure;
t = tiledlayout('flow')
for n = 1:nMic    % For each microphone signal
    % Load the signal

    
    % Comput the impulse response using the function extractirsweep
    [ir] =

    % Setting up time scale for computed ir
    t = 
    
    % Find the first impulse of the impulse response
    

    % Double check if we are finding back the correct direct path length
    directPathTimeOfArrival = 
    directPathLength(n) = 
    
    % Plot the estimated impulse response
    nexttile
    plot(t, ir);
    xlim([0 0.05]);
    xlabel('Time (sec)');
    title(['Mic: ', num2str(n)]);
    
end

%% SPEAKER TO MIC DISTANCE (DIRECT PATH LENGTH)
% Print on screen the estimated distance from the source
figure;
plot(1:nMic, directPathLength)
xlabel('Measurement'), ylabel('Distance highest peak')

fprintf(sprintf('Direct path length %f m\n', directPathLength));

%% MIC TO WALLS DISTANCE COMPUTATION (ROOM SIZE)
% Put here the difference between first reflection and direct sound time of
% arrivals (TOA)

% Inspecting the impulse responses determine the delay of the first
% reflection
delay = ; %[s]
% Compute the distance from the reflectors
distance = ;

% Print on screen the estimated distance from the walls
fprintf(sprintf('Average distance between first path and reflector %f m\n', ...
    distance));
