clc, clear all, close all;
%% Question2
L = 0.645;   % [m] Length of the string
fE2 = 82;   % [Hz]
mu = 5.78e-3;    % [kg m^-1] linear density
T_str = 4*(L^2)*mu*fE2^2;   % [N] tension
c = sqrt(T_str/mu);
beta = 1/5;

Fs = 44100; %sampling frequency
T_E2E1 = (2*L*(1-beta))/c;
T_E1R1 = (beta*L)/c;
T_S = (2*L)/c;

N_E = round(T_E2E1*Fs);
N_E1R1 = round(T_E1R1*Fs);
N_S = round(T_S*Fs);

%filter from Simulink
signalLen = 10;

%load simulation
open_system(['ques2.slx'], 'loadonly');

%perform simulation
simulation = sim(['ques2.slx'], signalLen);

input = simulation.input.data;
output = simulation.filter.data;

%bridge impedance from q1
% Parameters of the bridge 
k_p = 1.41e5; % [N/m] top plate stiffness
m_p = 0.128*0.385;    % [kg] top plate effective mass
A_p = 0.0375*0.385;   % [m^2] top plate area
R_p = 32; % [N s m^-1] top plate resistance

m_h = 8.04e-4;   % [kg]  air piston mass
A_h = 7.85e-3;   % [m^2] air piston area
R_h = 30;    % [Ns/m^5] sound hole resistance

V = 0.0172;  %[m^3] air cavity volume
rho = 1.2;   % [kg m^-3] air density
c = 343;    % [m/s] sound velocity in air

L_p = m_p/(A_p^2);  % plate inertance
L_h = m_h/(A_h^2);  % sound hole inertance
C_p = A_p^2/k_p;    % plate compliance
C_v = V/(rho*c^2);  % cavity compliance

%plot of the bridge impedance

f = 0:Fs/length(input):Fs-(1/length(input));  %[Hz] range of frequency in Hz
w = 2*pi*f; %[rad/s] range of angular frequency 

Z_p = 1i.*w*L_p + R_p + (1i.*w*C_p).^(-1); %[Kg/m^4s] plate impedance
Z_v = 1./(1i.*w*C_v); %[Ks/m^4 s] cavity impedance
Z_h = 1i.*w*L_h + R_h; %[Ks/m^4 s] hole impedance

Z_b = Z_p + ( Z_v .* Z_h ) ./ (Z_v + Z_h ); %[Ks/m^4 s] bridge impedance

%filter without bridge impedance
H_1 = (abs(fft(output) ./ fft(input)))';

%filter with bridge impedance
H = H_1.*Z_h;

semilogy(f, abs(H));
xlim([0 500])