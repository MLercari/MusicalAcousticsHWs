%% Es. 1 HL2
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

%resonance frequency
fd = omegad/(2*pi); 


%% Simulation

% Simulation Parameters
Fs = 100000;                  % Sampling Frequency
signalLen = 10;              % Simulation Duration

%load simulation
open_system(['Es1.slx'], 'loadonly');

%perform simulation
simulation = sim(['Es1.slx'], signalLen);

%Compute Frequency Response
input = simulation.input.data;
output = simulation.output.data;
f = 0:Fs/length(input):Fs-(1/length(input));    
H = abs(fft(output) ./ fft(input));

 %% compare with analytical frequency response
Z = ((1i*2*pi.*f).^2.*(M*C) + 1i*2*pi.*f.*(R*C) + 1)./(1i*C*2*pi.*f);
H_a = abs(Z.^(-1)); 
 
figure(1)
%plot of transfer function from simulated data

semilogx(f, 20*log10(H.*R));
pause(0.05);
hold on

xline(real(fd),'-',{'54522 Hz'});
%plot of tranfer function from analytical solution
semilogx(f, 20*log10(H_a.*R));
legend(["FRF from simulated data" , "Resonance frequency ", "FRF from analytical data" ],'Location','northwest')
xlabel("f [Hz]" ,'FontSize',12,'FontWeight','bold','Color','k','interpreter','latex')
ylabel("$|H(f) \times R|$ [$dB$]" ,'FontSize',12,'FontWeight','bold','Color','k','interpreter','latex')

xlim([0 50000])
grid on

