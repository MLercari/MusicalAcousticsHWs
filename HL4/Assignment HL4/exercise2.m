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

nMic = 24; % Number of microphones: assumed to be equal to the number of measurements

typeOfSignal = ["sweep" , "noise"]; % Noise or sweep
dir =  "Recordings";     % Recordings directory

c = 343.8; % [m]/[s]

%% Plot the signal autocorrelation

for n = 1:nMic
    
    % Load the signal
    [x, Fs] = audioread(strcat(dir,'/', typeOfSignal(2), '/', num2str(n) , '.wav' ));
    
    % Time length
    tlen = 4; %[s] duration of the signal
    tsamples = length(x);
    
    % Time axis
    t = linspace(0, tlen, tsamples)';
    
    % Auto correlation (Half of the autocorrelation)
    % xd = circshift(x,1);
    xc = xcorr(x, 'normalized');
    xc(1:tsamples-1) = [];

    % Plot the autocorrelations
    
    %{
    figure(n)
    plot(t, xc);
    title(['Mic: ', num2str(n)]);
    axis([0 0.05 -1 1]);    % Limit the axis
    xlabel('Time (sec)');
    %}
    
end

%% MIC TO REFLECTORS DISTANCE COMPUTATION

% Put here the difference between first reflection 

delay = zeros(nMic, 1); %every first reflection delay is stored here


%by visual inspection of the autocorrelations graphs we found that first
%reflection peak lies down around 3.5 ms, we find it by looking around a
%time interval [25 ms, 40 ms]


for jj = 1:nMic %compute the delay for every single measurement
    
    [x, Fs] = audioread(strcat(dir,'/', typeOfSignal(2), '/', num2str(jj) , '.wav' ));
    
    trunc_xc = xcorr(x, 'normalized');
    
    %remove half correlation
    trunc_xc(1:tsamples-1) = [];
    
    %after having found peaks position by visual inspection, we search
    %their excact values 

    leftoffset = 0.0025; %[s]
    rightoffset = 0.0040; %[s]
    trunc_xc(1:leftoffset*Fs) = [];  
    trunc_xc(rightoffset*Fs:end) = [];
    
    %search for the maxima
    [maxVal,maxIndex] = max(trunc_xc);
    
    %store the delay in a matrix
    delay(jj) = maxIndex/Fs; % [s] 
    
    %print them for inspection 
    time = linspace(leftoffset, length(trunc_xc)/Fs , length(trunc_xc));
    figure(jj)
    plot(time, trunc_xc)
    hold on
    xline(delay(jj)) 
    title(strcat("Mes n." , num2str(jj)))
    xlim([leftoffset rightoffset]);    % Limit the axis
    
end


distance = mean(delay)*c; %[m]

fprintf(sprintf('Average distance between mic and first reflection %f m\n', distance));
