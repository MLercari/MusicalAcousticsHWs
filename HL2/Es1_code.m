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

omega0 = 1/(sqrt(M*C));
% considering damping
omegad = omega0*sqrt(1 - (R/(2*M*omega0))^2); 
fd = omegad/(2*pi);


%%
% Simulation Parameters
Fs = 100000;                  % Sampling Frequency
signalLen = 10;              % Simulation Duration

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
    hold on
    xline(real(fd));