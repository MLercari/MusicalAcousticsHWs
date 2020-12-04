%% Es. 2 HL2
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

% current is measured at the end of each tree.
% Simulation Parameters
Fs = 100000;                  % Sampling Frequency
signalLen = 2;              % Simulation Duration

%N is the height and K are the number of leaves

indexes = [12, 21, 22, 23, 32, 33]; % NK
plotcolours = ["k", "m", "b", "g", "r", "c"];

for i = 1:6

%load simulation
open_system(["RT" + int2str(indexes(i)) + ".slx"], 'loadonly');

%perform simulation
simulation = sim(["RT" + int2str(indexes(i)) + ".slx"], signalLen);

%Compute Frequency Response
input = simulation.input.data;
output = simulation.output.data;
f = 0:Fs/length(input):Fs-(1/length(input));    
H = abs(fft(output) ./ fft(input));
    
figure(1)
subplot (3, 2, i)
plot(f, 20*log10(H.*R) , plotcolours(i) , "LineWidth", 1.2 );
pause(0.05);

xlabel("f [Hz]" ,'FontSize',12,'FontWeight','bold','Color','k','interpreter','latex')
ylabel("$|H(f)|\times R [dB]$ " ,'FontSize',12,'FontWeight','bold','Color','k','interpreter','latex')
xlim([0 50000])
ylim([ -60 0 ])
legend(int2str(indexes(i))); 
grid on

end

%legend(["NK 12", "NK 21", "NK 22", "NK 23", "NK 32", "NK 33"])
    
%% compare the simulated FRF with analytical FRF
 
% N = 2; %height
% K = 3; %leaves
%   
% Zs = ((1i*2*pi.*f).^2.*(L*C) + 1i*2*pi.*f.*(R*C) + 1)./(1i*C*2*pi.*f);
% 
%   
%   for n = 1:N
%       Z = 1i*2*pi*L.*f + R + 1./(1i*2*C*pi.*f + K./Zs);
%       Zs = Z;
%   end
%  
% H_a = abs(Z.^(-1));
%   
% figure(2)
% semilogx(f, 20*log10(H.*R));
% pause(0.05);
% hold on
% 
% %plot of tranfer function from analytical solution
% semilogx(f, 20*log10(H_a.*R));
% 
% legend(["FRF from simulated data 2x3" , "FRF from analytical data 2x3" ],'Location','northwest')
% xlabel("f [Hz]" ,'FontSize',12,'FontWeight','bold','Color','k','interpreter','latex')
% ylabel("$|H(f) \times R|$ [$dB$]" ,'FontSize',12,'FontWeight','bold','Color','k','interpreter','latex')
% 
% xlim([0 50000])
% grid on
    
