%% Es. 1 HL1
clar all;
close all;
clc;

%%
% Simulation Parameters
Fs = 2000;                  % Sampling Frequency
signalLen = 3;              % Simulation Duration

%load simulation
open_system(['Es1.slx'], 'loadonly');

%perform simulation
