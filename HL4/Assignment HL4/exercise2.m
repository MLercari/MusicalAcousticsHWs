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

nMic = 1; % Number of microphones: assumed to be equal to the number of measurements

typeOfSignal = ["sweep" , "noise"]; % Noise or sweep
duration = [10 , 4];
dir =  "Recordings";     % Recordings directory
Fs = 48000; %[Hz] Sampling frequency

c = 343.8; % [m]/[s]

%% Plot the signal autocorrelation

figure(1);
tiledlayout('flow');

%initialize the correlation vector
corr_sweep = zeros(nMic, Fs*10);
corr_noise = zeros(nMic, Fs*4);
corr = [{corr_sweep} ; {corr_noise}];

for m = 1:length(typeOfSignal)
    
    for n = 1:nMic     
    % Load the signal
    [x, Fs] = audioread(strcat(dir,'/', typeOfSignal(m), '/', num2str(n) , '.wav' ));
    x = x(1:Fs*duration(m));
    
    % Time length
    tlen = duration(m); %[s] duration of the signal
    tsamples = length(x);
    
    % Time axis
    t = linspace(0, tlen, tsamples)';
    
    % Auto correlation (Half of the autocorrelation)
    xc = xcorr(x, 'normalized');
    xc(1:tsamples-1) = [];

    % Plot the autocorrelations
    nexttile
    plot(t, xc);
    title(['Mic: ', num2str(n), typeOfSignal(m)]);
    %axis([0 0.05 -1 1]);    % Limit the axis
    xlim([ 0.0025 0.004])
    xlabel('Time (sec)');
  
    
    %store the autocorrelations
    current_corr = corr{m};
    current_corr(n , :) = xc;
    corr{m} = current_corr;
    end
end

%% MIC TO REFLECTORS DISTANCE COMPUTATION

% Put here the difference between first reflection 

delay = zeros(nMic, 1); %every first reflection delay is stored here

%by visual inspection of the autocorrelations graphs we found that first
%reflection peak lies around 3.5 ms (corresponding to the first reflection time),

figure(2);
tiledlayout('flow');

for ii = 1:length(typeOfSignal)
    
for jj = 1:nMic %compute the delay for every single measurement
    
    [x, Fs] = audioread(strcat(dir,'/', typeOfSignal(ii), '/', num2str(jj) , '.wav' ));
    x = x(1:Fs*duration(ii));
    
    trunc_xc = xcorr(x, 'normalized');
    
    % Time length
    tlen = duration(ii); %[s] duration of the signal
    tsamples = length(x);
    
    % Time axis
    t = linspace(0, tlen, tsamples)';
            
    %remove half correlation
    trunc_xc(1:tsamples-1) = [];
    
    %inspect the signal around 3.5 ms and find the peak position
    leftoffset = 0.0025; %[s]
    rightoffset = 0.0040; %[s]
    trunc_xc(1:leftoffset*Fs) = [];
    trunc_xc(rightoffset*Fs:end) = [];
    
    %time axis of the truncated vector
    trunc_t = t;
    trunc_t(1:leftoffset*Fs) = [];
    trunc_t(rightoffset*Fs:end) = [];
    
    %search for the maxima
    %[pks, locs] = findpeaks(trunc_xc, (leftoffset*Fs:rightoffset*Fs) , 'NPeaks', 1);
    
    %store the delay in a matrix
    %delay(jj) = locs/Fs; % [s] 
    
    %print them for inspection 
    nexttile
    plot(trunc_t, trunc_xc)
    title(strcat("Mes n." , num2str(jj) , typeOfSignal(ii)))
    xlim([leftoffset rightoffset]);    % Limit the axis
    %plot the peak 
    %{
    hold on
    y_stem = NaN(1,length(trunc_xc));
    y_stem(locs) = trunc_xc(locs);
    stem( time , y_stem)
    %plot(locs./Fs, pks, "LineStyle", "none", "Marker", "x", "Color", "m");
    

    hold off
    %}
    
end

end


distance = mean(delay)*c; %[m]

fprintf(sprintf('Average distance between mic and first reflection %f m\n', distance));
