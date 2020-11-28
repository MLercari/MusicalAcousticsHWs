%% Es. 1 HL1
clear all;
close all;
clc;

%physical parameters
V = 0.1; 
l = 0.1;
S = 100; 
rho = 1.2;
c = 343;

C = V/(rho*c^2); 
R = (rho*c)/S;
M = rho*l/S;

f0 = 1/(2*pi*sqrt(M*C));
% f0 is 5459 Hz, so sampling frequency must be at least 11000 Hz

%%
% Simulation Parameters
Fs = 20000;                  % Sampling Frequency
signalLen = 3;              % Simulation Duration

%load simulation
open_system(['Es1.slx'], 'loadonly');

%perform simulation
simulation = sim(['Es1.slx'], signalLen);

 % Compute Frequency Response
    input = simulation.input.data;
    output = simulation.output.data;
    f = 0:Fs/length(input):Fs-(1/length(input));    
    H = db(abs(fft(output) ./ fft(input)));
    
    plot(f, H);
    pause(0.05);