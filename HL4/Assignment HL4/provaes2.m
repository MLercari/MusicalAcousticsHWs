%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Homework 4 exercise 2
% Reflection inspection
% By a visual inspection of the recorded signals we evaluate the first
% reflection time from which the distance of the micrphones from the
% reflectors can be inferred.
%
% Musical Acoustic Course
% Mirco Pezzoli
% 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc 

%% First reflection estimation using autocorrelation
% Collecting autocorrelations and plotting them in order to analyze
% reflections. Try different takes and input signals

addpath('Functions');

nMic = 24; % Number of microphones
typeOfSignal = ['noise'; 'sweep']; % Noise or sweep
dir = './Recordings';      % Recordings directory

speed_of_sound = 343.8; % [m]/[s]

%% Plot the signal autocorrelation
figure;
t = tiledlayout('flow');
corr = zeros(nMic, 0.02*48e3);

for n = 1:nMic
    % Load the signal
    fullFileName = fullfile(dir,typeOfSignal(1,:),[num2str(n) '.wav']); 
    [x, fs] = audioread(fullFileName);
    
    % Time length
    dur = length(x)/fs;
      
    % Time axis
    t = 0:1/fs:dur-1/fs;
    
    % Auto correlation
    
    [xc, lags] = xcorr(x, 'normalized');
    
    %Half of autocorrelation
    
    xc = xc(lags >= 0);
    
    %Save a portion of autocorrelation to compute First reflection
    select = round(fs*0.002):round(fs*0.0044);
    corr(n,select) = xc(select);
    
    % Plot the autocorrelation
    nexttile
    plot(t, xc);
    title(['Mic: ', num2str(n)]);
    axis([0 0.02 -1 1]);    % Limit the axis
    xlabel('Time (sec)');
  
end

%% MIC TO REFLECTORS DISTANCE COMPUTATION
% Put here the difference between first reflection 

delay = zeros(1, 24);

for n = 1:nMic
    %select the portion of autocorrelation that contains the first
    %reflection peak
    select = corr(n,round(fs*0.002):round(fs*0.044));
    %find the maximum index of the maximum
    [~,I] = max(select);
    %compute the delay (remembering to offset by 2 ms)
    delay(n) = I/fs + 0.002;%[s]
end

delay_avg = sum(delay)/nMic;% [s]
distance = sum(speed_of_sound*delay)/24; %[m]

fprintf(sprintf('Average distance between mic and first reflection %f m\n', distance));
