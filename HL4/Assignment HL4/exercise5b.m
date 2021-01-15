%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Homework 4 exercise 5b
% Radiance pattern estimation from impulse response
% We compute the impulse response using the sine sweep source signal
% from which the radiance pattern is estimated windowing the signal befor
% the occurence of the first reflection
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
addpath("Recordings");

fs = 48000;         % Sampling frequency
nfft = fs;          % Number of fft points
duration = 10;       % [s] duration of sweep signal
nMic = 24;          % Number of microphones
R = 2.67;            % Distance between source and microphones

c = 343.8; % [m]/[s]

typeOfSignal = "sweep"; % Sweep
dir = strcat("Recordings/", typeOfSignal, "/"); % File directory

%% RADIANCE COMPUTATION THROUGH THE USE OF THE IMPULSE RESPONSE


% First reflections attenuation. We consider a small opart of the ir
% windowing the impulse response
ns = 50;        % We need a small window
%create a hanning window of 2*ns +1 length
hann_win =  hann(2*ns +1);     % Window



TOA_directSignal = 0.00777;     % TOA
TOA_firstReflection =  0.010983; % First reflection TOA



% This is the input signal. The following function computes 1/FFT(sweep). 
% We will thus need to MULTIPLY this signal for FFT(output) in order to 
% apply deconvolution
% Use the provided synthSweep function with frequency interval between 50Hz
% and 22kHz
[x invsweepfft sweepRate] = synthSweep(duration, fs, 50, 22e3);


figure;
tiledlayout('flow');
sig = zeros(fs*duration/2, nMic);

for n = 1:nMic                    % For each microphone signal
    % Load the signal
    y = audioread(strcat(dir, num2str(n), ".wav"));
    y = y(1:fs*duration);
    
    % Compute the impulse response using the function extractirsweep
    [ir] = extractirsweep(y, invsweepfft);
    
    % Setting up time scale for computed ir
    t = (0:1/fs:length(ir)/fs);
    t = t(1:end -1);
    
    % Windowing the ir. In order to take into account windows of arbitrary
    % length, the ir is circularly shifted until the max peak is at the
    % middle of the array. This prevents the left boundary of the windows
    % to be out of the vector. See function circshift.
       
    % Use the TOA information in order to find the correct poisition of the
    % first path in the shifted signal.
    
    % Impulse response windowing
    % we want to maintain time information in the windowed impulse response.
    % We allocate one array of the same length of the original ir, and then
    % add samples from the windowed signal. Finally, we invert the
    % previously applied circShift to return back to the correct impulse
    % response timing.
    
    shiftAmount = round(length(ir)/2 - TOA_directSignal*fs);
    irCircular = circshift(ir, shiftAmount);  
    w = zeros(length(ir),1);
    w(length(ir)/2-ns:length(ir)/2+ns) = hann_win;
    irCircular_w =  irCircular.*w; % Window the ir
    ir_w = circshift(irCircular_w , - shiftAmount); % Shift ir back
    
    
   % Plot the estimated impulse response
    nexttile
    hold on
    plot(t, ir, "Color", "m")

    % Plot the TOA and TOA First over the IR
    x_stem = t;
    y_stem1 = NaN(length(t),1);
    y_stem1(round(TOA_directSignal*fs)) = ir(round(TOA_directSignal*fs));
    stem(x_stem, y_stem1 );
    y_stem2 = NaN(length(t),1);
    y_stem2(round(TOA_firstReflection*fs)) = ir(round(TOA_firstReflection*fs));
    stem(x_stem, y_stem2 );
    
    % Plot the windowed IR over with thicker line
    plot(t, ir_w, 'LineWidth' , 0.5)
    
    %plot the window
    plot(t, circshift(w, -shiftAmount).*0.008);
    
    hold off
    
    xlim([0 0.015]);
    xlabel('Time (sec)');
    title(['Mic: ', num2str(n)]);
    
    sig(:,n) = ir_w;
end

%% Radiance estimation

SIG = fft(sig, nfft); % FFT of the windowed signal

f = 0:fs; 
f = f(1:end -1);

G = exp(-1i*(f./c)*R)./(4*pi*R); %green function

rad_patt = zeros(fs, nMic); % Compute the radiance pattern magnitude

for i = 1:nMic
    
    rad_patt(:,i) = SIG(:,i)./G'; % Compute the radiance pattern magnitude

end

%% Radiance pattern plot
% Frequencies must be centered at a frequency bins
ctr_freqs = [150 250 500 750 1000 1500 2000 5000 7500 10000];

% Microphone angle wrt the speaker frontal face
angs = zeros(1,nMic);

for ii = 2:nMic
    angs(1,ii) = angs(1,ii-1) + 15;
end

% Plot the estimated radiance pattern using the provided function
% radianceplot

radianceplot(ctr_freqs, rad_patt, angs, 'sweep IR: ');

