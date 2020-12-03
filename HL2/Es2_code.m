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
open_system(['subsToComp.slx'], 'loadonly');

%perform simulation
simulation = sim(['subsToComp.slx'], signalLen);

 % Compute Frequency Response
    input = simulation.input.data;
    output = simulation.output.data;
    f = 0:Fs/length(input):Fs-(1/length(input));    
    H = abs(fft(output) ./ fft(input));
    


    
 %%
 % compare with analytical frequency response
N = 2; %height
K = 3; %leaves
  
 Zs = ((1i*2*pi.*f).^2.*(L*C) + 1i*2*pi.*f.*(R*C) + 1)./(1i*C*2*pi.*f);

  
  for n = 1:N
      Z = 1i*2*pi*L.*f + R + 1./(1i*2*C*pi.*f + K./Zs);
      Zs = Z;
  end
 

%   for n = 1:N-1
%       Z = ((1i*2*pi.*f).^2.*(L*C) + 1i*2*pi.*f.*(C*R) + 1)./(1i*2*C*pi.*f + K./Zs0);
%       Zs = Z;
%       Z = 1i*2*pi*L.*f + R + 1./(1i*2*C*pi.*f + 1./Zs + 1./Zs0);
%   end

  H_a = abs(Z.^(-1));
  
   figure(1)
    plot(f, H);
    pause(0.05);
    hold on
    plot(f, H_a);
    legend(["simulation" , "analytical"]);
    grid on
    xlim([0 50000])
    

