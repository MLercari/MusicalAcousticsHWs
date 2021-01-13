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

nMic = 1;% Number of microphones
typeOfSignal = ["sweep" , "noise"]; % Noise or sweep
dir =  "Recordings";     % Recordings directory

c = 343.8; % [m]/[s]

%% Plot the signal autocorrelation
figure;
%t = tiledlayout('flow');

for n = 1:nMic
    
    % Load the signal
    [x, Fs] = audioread(strcat(dir,'/', typeOfSignal(2), '/', '1.wav' ));
    
    % Time length
    tlen = 4; %[s] duration of the signal
    tsamples = length(x);
    
    % Time axis
    t = linspace(0, tlen, tsamples)';
    
    % Auto correlation (Half of the autocorrelation)
    xd = circshift(x,1);
    xc = xcorr(x, 'normalized');
    xc(1:tsamples-1) = [];

    % Plot the autocorrelation
    nexttile
    plot(t, xc);
    title(['Mic: ', num2str(n)]);
    axis([0 0.05 -0.04 0.08]);    % Limit the axis
    xlabel('Time (sec)');
    
end

%% MIC TO REFLECTORS DISTANCE COMPUTATION

% Put here the difference between first reflection 

nMes = 24; %number of measurements

delay = zeros(nMes, 1); %every first reflection delay is stored here
ampPeak = zeros(nMes, 1); %the relative first peak amplitude

for jj = 1:nMes %compute the delay for every single measurement
    
    [x, Fs] = audioread(strcat(dir,'/', typeOfSignal(2), '/', num2str(jj) , '.wav' ));
    
    trunc_xc = xcorr(x, 'normalized');
    
    %remove half correlation
    trunc_xc(1:tsamples-1) = [];
    
    %remove initial spike and signal after 10 milli seconds
    offset = 0.001;
    trunc_xc(1:offset*Fs) = [];  
    trunc_xc(0.01*Fs:end) = [];
    
    %search for the maxima
    [maxVal,maxIndex] = max(trunc_xc);
    
    %store the delay in a matrix
    delay(jj) = maxIndex/Fs; % [s] 
    ampPeak(jj) = maxVal;
    
    %print them for inspection 
    time = linspace(0, length(trunc_xc)/Fs, length(trunc_xc));
    figure(jj)
    plot(time, trunc_xc)
    hold on
    xline(delay(jj) + offset) 
    title(strcat("Mes n." , num2str(jj)))
    axis([offset 0.01 -0.05 0.05]);    % Limit the axis
    
end

distance = mean(delay)*c; %[m]

fprintf(sprintf('Average distance between mic and first reflection %f m\n', distance));
