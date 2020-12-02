%% Es. 2 HL1
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
L = rho*l/S;

%% Simulation

% Simulation Parameters
Fs = 100000;                  % Sampling Frequency
signalLen = 10;              % Simulation Duration

%load simulation
open_system(['Es2.slx'], 'loadonly');

%perform simulation
simulation = sim(['Es2.slx'], signalLen);

 % Compute Frequency Response
    input = simulation.input.data;
    output = simulation.output.data;
    f = 0:Fs/length(input):Fs-(1/length(input));    
    H = db(abs(fft(output) ./ fft(input)));
    
    plot(f, H);
    pause(0.05);
    hold on
    xline(real(fd));
    
 %%
 % compare with analytical frequency response

