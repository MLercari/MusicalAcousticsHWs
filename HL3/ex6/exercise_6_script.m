%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 6
% Complete model of a guitar
% This script call the simulink file that implements the simple electrical
% equivalent of a guitar. The simulation is executed and the resulting
% current is resampled with the desired sampling frequency. Finally the
% sound is plotted in time and frequency domain and saved on disk. 
%
% Musical Acoustic Course
% Mirco Pezzoli
% 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;

addpath('Functions')
simulink_folder = './';       % Simulink projects folder
addpath(simulink_folder);

%% Setup
fs = 44100;                         % Sampling frequency
signalLen = 3;                      % Signal length
t = [0:1/fs:signalLen-1/fs];        % Time axis

fileName = 'complete_guitar.wav';     % Audio file path

%% Simulation
% run the simulink simulation using the command sim (see doc sim).
out = sim('exercise_6');
%dt = diff(out.I.time);
%force = timeseries((diff(out.I.data)).*dt, out.I.time(1:end-1));

% The variable I contains non constant time intervals between samples.
% Resample the data using resample function in order to obtain an equally
% sampled signal I1
I1 = resample(out.I, t);

% Plot the resampled signal in time
figure(1);
plot(I1.time, I1.data);


% Normalize the signal
soundWave = I1.data;
soundWave = diff(soundWave)*fs;
soundWave = soundWave./max(abs(soundWave));


%% Plot and play
% Plot the signal frequency content as magnitude and phase

[S, magS, angleS, f, df] = myFFT(soundWave, fs);

h = figure(2);
plotFFT_linearFreqScale(magS, angleS, f, df, fs, 2000, h);


figure(3);
t(end) = [];
plot(t, soundWave);
title('Signal in time'), xlabel('Time [s]')

sound(soundWave, fs);                           % Play the sound

disp('Save file on disk...')                    % Save on disk
audiowrite(fileName, soundWave, fs);
