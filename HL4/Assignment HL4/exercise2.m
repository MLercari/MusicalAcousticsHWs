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
duration = [10 , 3];
dir =  "Recordings";     % Recordings directory
Fs = 48000; %[Hz] Sampling frequency

c = 343.8; % [m]/[s]

%% Plot the signal autocorrelation

figure(1);
tiledlayout('flow');

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
  
    end
end

%% MIC TO REFLECTORS DISTANCE COMPUTATION

% Put here the difference between first reflection 

delay = zeros(nMic, 2); %every first reflection delay is stored here

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

    trunc_xc = trunc_xc(leftoffset*Fs:rightoffset*Fs);   
      
    %time axis of the truncated vector
    trunc_t = t;

    trunc_t = trunc_t(leftoffset*Fs:rightoffset*Fs); 
    %search for the maxima
    [pks, locs] = findpeaks(trunc_xc, trunc_t , 'NPeaks', 1 , 'SortStr' , ...
        'descend' , 'MinPeakHeight' , -1);
    
    %sometimes the peak detector doesn't work and it returns an empty
    %vector: to avoid this we set 0.0035 
    if (isempty(locs) ==1)
        locs = 0.0035;
    end
    
    %store the delay in a matrix
    delay(jj , ii) = locs; % [s] 
    
    %print them for inspection 
    nexttile
    plot(trunc_t, trunc_xc)
    title(strcat("Mes n." , num2str(jj) , typeOfSignal(ii)))
    xlim([leftoffset rightoffset]);    % Limit the axis
    %plot the peak 
    hold on
    y_stem = NaN(1,length(trunc_t));
    y_stem( round((locs - leftoffset)*Fs)) = trunc_xc( round((locs - leftoffset)*Fs));
    stem( trunc_t , y_stem)
    
    hold off

    
end

end


distance = mean(delay)*c; %[m]
distance = mean(distance'); 

mean_delay = mean(delay);
variance_delay = var(delay);

fprintf(sprintf('Average distance between mic and first reflection %f m\n', distance));
fprintf(sprintf('with a mean time of first reflection %f s\n', mean_delay));
fprintf(sprintf('and relative variance %f s\n' , variance_delay ));
